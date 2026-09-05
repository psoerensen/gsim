#include "standalone_vcf.h"

#include "standalone_metadata_error.h"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace gsim::native::metadata {

class VcfLineInput {
 public:
  virtual ~VcfLineInput() = default;
  virtual bool next(std::string& line) = 0;
  virtual void finish() = 0;
  [[nodiscard]] virtual std::uint64_t fixed_buffer_bytes() const noexcept = 0;
};

namespace {

namespace fs = std::filesystem;
constexpr std::size_t kCompressedBufferBytes = 65536u;

[[noreturn]] void invalid(const std::string& message) {
  throw Error(StatusCode::invalid_argument, "VCF: " + message);
}

[[noreturn]] void io_failure(const std::string& message) {
  throw Error(StatusCode::internal_error, "VCF: " + message);
}

std::vector<std::string> fields(const std::string& value, char separator) {
  std::vector<std::string> result;
  std::size_t begin = 0;
  while (true) {
    const auto end = value.find(separator, begin);
    result.push_back(value.substr(begin, end - begin));
    if (end == std::string::npos) break;
    begin = end + 1u;
  }
  return result;
}

std::uint64_t position(const std::string& value, std::uint64_t line) {
  if (value.empty() || value[0] == '-' || value[0] == '+') {
    invalid("invalid POS at line " + std::to_string(line));
  }
  errno = 0;
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
  if (errno == ERANGE || end == value.c_str() || *end != '\0' || parsed == 0) {
    invalid("invalid POS at line " + std::to_string(line));
  }
  return static_cast<std::uint64_t>(parsed);
}

std::size_t gt_index(const std::string& format, std::uint64_t line,
                     std::size_t& format_count) {
  const auto values = fields(format, ':');
  format_count = values.size();
  std::size_t found = values.size();
  std::unordered_set<std::string> names;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (values[i].empty() || !names.insert(values[i]).second) {
      invalid("FORMAT fields must be nonempty and unique at line " +
              std::to_string(line));
    }
    if (values[i] == "GT") found = i;
  }
  if (found == values.size()) {
    invalid("FORMAT does not contain GT at line " + std::to_string(line));
  }
  return found;
}

enum class Unsupported {
  none,
  indel,
  multiallelic,
  symbolic,
  other_allele,
  missing,
  unphased,
  non_diploid
};

const char* unsupported_name(Unsupported value) {
  switch (value) {
    case Unsupported::indel: return "indel";
    case Unsupported::multiallelic: return "multiallelic ALT";
    case Unsupported::symbolic: return "symbolic or breakend allele";
    case Unsupported::other_allele: return "unsupported or identical allele";
    case Unsupported::missing: return "missing selected-sample GT";
    case Unsupported::unphased: return "unphased selected-sample GT";
    case Unsupported::non_diploid: return "non-diploid selected-sample GT";
    case Unsupported::none: break;
  }
  return "supported";
}

bool is_base(const std::string& value) {
  return value.size() == 1u &&
         (value[0] == 'A' || value[0] == 'C' ||
          value[0] == 'G' || value[0] == 'T');
}

Unsupported allele_reason(const std::string& reference,
                          const std::string& alternate) {
  if (alternate.find(',') != std::string::npos) return Unsupported::multiallelic;
  const auto symbolic = [](const std::string& value) {
    return value.find('<') != std::string::npos ||
           value.find('>') != std::string::npos ||
           value.find('[') != std::string::npos ||
           value.find(']') != std::string::npos || value == "*";
  };
  if (symbolic(reference) || symbolic(alternate)) return Unsupported::symbolic;
  if (reference.size() != 1u || alternate.size() != 1u) return Unsupported::indel;
  if (!is_base(reference) || !is_base(alternate) || reference == alternate) {
    return Unsupported::other_allele;
  }
  return Unsupported::none;
}

bool numeric_token(const std::string& token) {
  return !token.empty() &&
         std::all_of(token.begin(), token.end(),
                     [](unsigned char value) { return std::isdigit(value) != 0; });
}

Unsupported classify_gt(const std::string& value, std::uint8_t* left,
                        std::uint8_t* right, std::uint64_t line) {
  if (value.empty()) invalid("empty GT at line " + std::to_string(line));
  for (const unsigned char value_byte : value) {
    if (!(std::isdigit(value_byte) || value_byte == '.' ||
          value_byte == '|' || value_byte == '/')) {
      invalid("malformed GT syntax at line " + std::to_string(line));
    }
  }
  const bool has_pipe = value.find('|') != std::string::npos;
  const bool has_slash = value.find('/') != std::string::npos;
  if (has_pipe && has_slash) {
    invalid("GT mixes phased and unphased separators at line " +
            std::to_string(line));
  }
  const char separator = has_pipe ? '|' : (has_slash ? '/' : '\0');
  if (separator == '\0') {
    if (value == ".") return Unsupported::missing;
    if (!numeric_token(value)) invalid("malformed GT syntax at line " + std::to_string(line));
    return Unsupported::non_diploid;
  }
  const auto alleles = fields(value, separator);
  for (const auto& allele : alleles) {
    if (allele.empty()) invalid("malformed GT syntax at line " + std::to_string(line));
    if (allele == ".") continue;
    if (!numeric_token(allele)) invalid("malformed GT syntax at line " + std::to_string(line));
  }
  if (std::any_of(alleles.begin(), alleles.end(),
                  [](const std::string& allele) { return allele == "."; })) {
    return Unsupported::missing;
  }
  if (alleles.size() != 2u) return Unsupported::non_diploid;
  if (separator == '/') return Unsupported::unphased;
  if ((alleles[0] != "0" && alleles[0] != "1") ||
      (alleles[1] != "0" && alleles[1] != "1")) {
    invalid("GT allele index is invalid for a biallelic record at line " +
            std::to_string(line));
  }
  if (left != nullptr) {
    *left = static_cast<std::uint8_t>(alleles[0][0] - '0');
    *right = static_cast<std::uint8_t>(alleles[1][0] - '0');
  }
  return Unsupported::none;
}

struct VcfRecord {
  std::array<std::string, 9> fixed;
  std::size_t sample_count{0};
  std::size_t sample_offset{0};

  [[nodiscard]] const std::string& operator[](std::size_t index) const {
    return fixed[index];
  }
};

VcfRecord record_fields(std::string& line, std::uint64_t number,
                        std::size_t sample_count) {
  if (!line.empty() && line.back() == '\r') line.pop_back();
  VcfRecord result;
  std::size_t begin = 0;
  std::size_t field = 0;
  while (true) {
    const auto end = line.find('\t', begin);
    if (field < result.fixed.size()) {
      result.fixed[field] = line.substr(begin, end - begin);
    }
    if (field == 8u) result.sample_offset = end + 1u;
    ++field;
    if (end == std::string::npos) break;
    begin = end + 1u;
  }
  if (field != 9u + sample_count) {
    invalid("sample-field count does not match #CHROM header at line " +
            std::to_string(number));
  }
  result.sample_count = sample_count;
  return result;
}

std::vector<std::string> selected_fields(
    const std::string& line, const VcfRecord& record,
    const std::vector<std::size_t>& selected_columns, std::uint64_t number) {
  std::vector<std::pair<std::size_t, std::size_t>> order;
  order.reserve(selected_columns.size());
  for (std::size_t output = 0; output < selected_columns.size(); ++output) {
    order.emplace_back(selected_columns[output], output);
  }
  std::sort(order.begin(), order.end());
  std::vector<std::string> result(selected_columns.size());
  std::size_t begin = record.sample_offset;
  std::size_t sample = 0;
  std::size_t wanted = 0;
  while (wanted < order.size()) {
    const auto end = line.find('\t', begin);
    if (sample == order[wanted].first) {
      result[order[wanted].second] = line.substr(begin, end - begin);
      ++wanted;
    }
    if (end == std::string::npos) break;
    begin = end + 1u;
    ++sample;
  }
  if (wanted != order.size()) {
    invalid("selected sample field is absent at line " + std::to_string(number));
  }
  return result;
}

std::string selected_gt(const std::string& sample_field, std::size_t gt,
                        std::size_t format_count, std::uint64_t line) {
  const auto values = fields(sample_field, ':');
  if (values.size() != format_count) {
    invalid("selected sample FORMAT arity mismatch at line " +
            std::to_string(line));
  }
  return values[gt];
}

Unsupported genotype_reason(const std::vector<std::string>& sample_fields,
                            std::size_t gt, std::size_t format_count,
                            std::uint64_t line) {
  bool missing = false;
  bool unphased = false;
  bool non_diploid = false;
  for (const auto& sample_field : sample_fields) {
    const Unsupported value = classify_gt(
        selected_gt(sample_field, gt, format_count, line), nullptr, nullptr, line);
    missing = missing || value == Unsupported::missing;
    unphased = unphased || value == Unsupported::unphased;
    non_diploid = non_diploid || value == Unsupported::non_diploid;
  }
  if (missing) return Unsupported::missing;
  if (unphased) return Unsupported::unphased;
  if (non_diploid) return Unsupported::non_diploid;
  return Unsupported::none;
}

void count_unsupported(VcfImportReport& report, Unsupported value) {
  switch (value) {
    case Unsupported::indel: ++report.indels; break;
    case Unsupported::multiallelic: ++report.multiallelic_records; break;
    case Unsupported::symbolic: ++report.symbolic_or_breakend_alleles; break;
    case Unsupported::other_allele: ++report.other_unsupported_alleles; break;
    case Unsupported::missing: ++report.missing_gt; break;
    case Unsupported::unphased: ++report.unphased_gt; break;
    case Unsupported::non_diploid: ++report.non_diploid_gt; break;
    case Unsupported::none: break;
  }
}

class PlainInput final : public VcfLineInput {
 public:
  explicit PlainInput(const std::string& path) {
    stream_.open(fs::u8path(path), std::ios::binary);
    if (!stream_.is_open()) io_failure("could not open input: " + path);
  }

  bool next(std::string& line) override {
    if (std::getline(stream_, line)) return true;
    if (stream_.bad()) io_failure("error while reading plain VCF");
    return false;
  }

  void finish() override {
    if (stream_.bad()) io_failure("error while reading plain VCF");
    stream_.close();
  }

  [[nodiscard]] std::uint64_t fixed_buffer_bytes() const noexcept override {
    return 0u;
  }

 private:
  std::ifstream stream_;
};

class GzipInput final : public VcfLineInput {
 public:
  explicit GzipInput(const std::string& path) : path_(path) {
#ifdef _WIN32
    file_ = gzopen_w(fs::u8path(path).c_str(), "rb");
#else
    file_ = gzopen(path.c_str(), "rb");
#endif
    if (file_ == nullptr) io_failure("could not open compressed input: " + path);
  }

  ~GzipInput() override {
    if (file_ != nullptr) (void)gzclose(file_);
  }

  bool next(std::string& line) override {
    line.clear();
    while (true) {
      for (std::size_t i = position_; i < available_; ++i) {
        if (buffer_[i] == '\n') {
          line.append(buffer_.data() + position_, i - position_);
          position_ = i + 1u;
          return true;
        }
      }
      if (available_ > position_) {
        line.append(buffer_.data() + position_, available_ - position_);
      }
      position_ = 0u;
      available_ = 0u;
      const int count = gzread(file_, buffer_.data(),
                               static_cast<unsigned int>(buffer_.size()));
      if (count < 0) compressed_error();
      if (count == 0) {
        int status = Z_OK;
        const char* detail = gzerror(file_, &status);
        if (!gzeof(file_) || (status != Z_OK && status != Z_STREAM_END)) {
          io_failure("compressed input is truncated or corrupt: " + path_ +
                     (detail == nullptr ? std::string() : ": " + std::string(detail)));
        }
        return !line.empty();
      }
      available_ = static_cast<std::size_t>(count);
    }
  }

  void finish() override {
    if (file_ == nullptr) return;
    const int status = gzclose(file_);
    file_ = nullptr;
    if (status != Z_OK) {
      io_failure("compressed input failed checksum/finalization validation: " + path_);
    }
  }

  [[nodiscard]] std::uint64_t fixed_buffer_bytes() const noexcept override {
    return static_cast<std::uint64_t>(buffer_.size());
  }

 private:
  [[noreturn]] void compressed_error() const {
    int status = Z_OK;
    const char* detail = gzerror(file_, &status);
    io_failure("compressed input decompression failed: " + path_ +
               (detail == nullptr ? std::string() : ": " + std::string(detail)));
  }

  std::string path_;
  gzFile file_{nullptr};
  std::array<char, kCompressedBufferBytes> buffer_{};
  std::size_t position_{0};
  std::size_t available_{0};
};

std::string detect_input_type(const std::string& path) {
  std::ifstream stream(fs::u8path(path), std::ios::binary);
  if (!stream.is_open()) io_failure("could not open input: " + path);
  std::array<unsigned char, 12> prefix{};
  stream.read(reinterpret_cast<char*>(prefix.data()),
              static_cast<std::streamsize>(prefix.size()));
  const std::streamsize count = stream.gcount();
  if (count < 2 || prefix[0] != 0x1fu || prefix[1] != 0x8bu) return "plain VCF";
  if (count < 12 || (prefix[3] & 0x04u) == 0u) return "gzip VCF";
  const std::uint16_t extra_length = static_cast<std::uint16_t>(prefix[10]) |
      (static_cast<std::uint16_t>(prefix[11]) << 8u);
  std::vector<unsigned char> extra(extra_length);
  if (extra_length != 0u) {
    stream.read(reinterpret_cast<char*>(extra.data()),
                static_cast<std::streamsize>(extra.size()));
    if (stream.gcount() != static_cast<std::streamsize>(extra.size())) {
      return "gzip VCF";
    }
  }
  std::size_t cursor = 0;
  while (cursor + 4u <= extra.size()) {
    const std::uint16_t length = static_cast<std::uint16_t>(extra[cursor + 2u]) |
        (static_cast<std::uint16_t>(extra[cursor + 3u]) << 8u);
    cursor += 4u;
    if (cursor + length > extra.size()) return "gzip VCF";
    if (extra[cursor - 4u] == 'B' && extra[cursor - 3u] == 'C' && length == 2u) {
      return "BGZF VCF";
    }
    cursor += length;
  }
  return "gzip VCF";
}

std::unique_ptr<VcfLineInput> open_input(const std::string& path,
                                         const std::string& input_type) {
  if (input_type == "plain VCF") return std::make_unique<PlainInput>(path);
  return std::make_unique<GzipInput>(path);
}

void unsupported_record(bool skip, VcfImportReport& report, Unsupported reason,
                        const VcfRecord& record,
                        std::uint64_t pos) {
  count_unsupported(report, reason);
  if (!skip) {
    invalid("unsupported record " + record[0] + ":" + std::to_string(pos) +
            " ID '" + record[2] + "': " + unsupported_name(reason));
  }
}

}  // namespace

PhasedVcfReader::PhasedVcfReader(
    const std::string& path, const std::vector<std::string>& selected_samples,
    const std::string& selected_chromosome, bool has_region,
    std::uint64_t region_start, std::uint64_t region_end,
    bool skip_unsupported)
    : path_(path) {
  if (path.empty()) invalid("path must not be empty");
  if (has_region && (selected_chromosome.empty() || region_start == 0u ||
                     region_start > region_end)) {
    invalid("region requires one chromosome and positive inclusive start <= end");
  }
  report_.input_type = detect_input_type(path_);
  auto input = open_input(path_, report_.input_type);
  std::unordered_set<std::string> header_sample_ids;
  std::unordered_set<std::string> retained_variant_ids;
  std::unordered_set<std::string> closed_chromosomes;
  std::string current_chromosome;
  std::uint64_t previous_position = 0;
  std::uint64_t line_number = 0;
  std::uint64_t maximum_line_bytes = 0;
  bool header = false;
  bool fileformat = false;
  bool saw_selected_chromosome = selected_chromosome.empty();
  std::vector<std::string> header_samples;
  std::string line;
  while (input->next(line)) {
    ++line_number;
    maximum_line_bytes = std::max(maximum_line_bytes,
                                  static_cast<std::uint64_t>(line.size() + 1u));
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!header) {
      if (line.rfind("##", 0) == 0) {
        if (line.rfind("##fileformat=VCFv", 0) == 0) fileformat = true;
        continue;
      }
      if (line.rfind("#CHROM\t", 0) != 0) {
        invalid("expected #CHROM header at line " + std::to_string(line_number));
      }
      const auto heading = fields(line, '\t');
      if (heading.size() < 10u || heading[0] != "#CHROM" ||
          heading[1] != "POS" || heading[2] != "ID" ||
          heading[3] != "REF" || heading[4] != "ALT" ||
          heading[8] != "FORMAT") {
        invalid("malformed #CHROM header");
      }
      if (!fileformat) invalid("missing ##fileformat=VCFv declaration");
      header_samples.assign(heading.begin() + 9, heading.end());
      for (std::size_t i = 0; i < header_samples.size(); ++i) {
        const auto& sample = header_samples[i];
        if (sample.empty() || sample == "." ||
            !header_sample_ids.insert(sample).second) {
          invalid("sample IDs must be unique, nonempty, and not '.'");
        }
      }
      report_.vcf_sample_count = static_cast<std::uint64_t>(header_samples.size());
      std::unordered_map<std::string, std::size_t> header_index;
      for (std::size_t i = 0; i < header_samples.size(); ++i) {
        header_index.emplace(header_samples[i], i);
      }
      if (selected_samples.empty()) {
        samples_ = header_samples;
        selected_columns_.resize(header_samples.size());
        for (std::size_t i = 0; i < header_samples.size(); ++i) selected_columns_[i] = i;
      } else {
        std::unordered_set<std::string> requested;
        for (const auto& sample : selected_samples) {
          if (sample.empty() || !requested.insert(sample).second) {
            invalid("selected sample IDs must be unique and nonempty");
          }
          const auto found = header_index.find(sample);
          if (found == header_index.end()) {
            invalid("selected sample '" + sample + "' is absent from the VCF header");
          }
          samples_.push_back(sample);
          selected_columns_.push_back(found->second);
        }
      }
      report_.selected_sample_count = static_cast<std::uint64_t>(samples_.size());
      header = true;
      continue;
    }
    if (line.empty() || line[0] == '#') {
      invalid("blank or header line after #CHROM at line " +
              std::to_string(line_number));
    }
    ++report_.total_records_scanned;
    const auto record = record_fields(line, line_number, header_samples.size());
    if (record[0].empty()) invalid("empty chromosome at line " + std::to_string(line_number));
    const std::uint64_t pos = position(record[1], line_number);

    if (record[0] != current_chromosome) {
      if (!current_chromosome.empty()) closed_chromosomes.insert(current_chromosome);
      if (closed_chromosomes.count(record[0]) != 0u) {
        invalid("chromosome '" + record[0] + "' occurs in disjoint blocks");
      }
      current_chromosome = record[0];
      previous_position = 0;
    }
    if (pos < previous_position) {
      invalid("POS decreases within chromosome '" + current_chromosome + "'");
    }
    previous_position = pos;

    if (!selected_chromosome.empty() && record[0] != selected_chromosome) {
      ++report_.outside_selected_chromosome;
      continue;
    }
    saw_selected_chromosome = true;
    if (has_region && (pos < region_start || pos > region_end)) {
      ++report_.outside_selected_region;
      continue;
    }
    if (record[2].empty()) invalid("empty variant ID at line " + std::to_string(line_number));
    std::string id = record[2];
    const bool generated = id == ".";
    if (generated) id = record[0] + ":" + record[1] + ":" + record[3] + ":" + record[4];
    if (!retained_variant_ids.insert(id).second) {
      ++report_.duplicate_final_ids;
      invalid("duplicate final variant ID '" + id + "'");
    }

    std::size_t format_count = 0;
    const std::size_t gt = gt_index(record[8], line_number, format_count);
    const auto selected_sample_fields =
        selected_fields(line, record, selected_columns_, line_number);
    for (const auto& sample_field : selected_sample_fields) {
      (void)selected_gt(sample_field, gt, format_count, line_number);
    }

    const Unsupported allele_status = allele_reason(record[3], record[4]);
    if (allele_status != Unsupported::none) {
      unsupported_record(skip_unsupported, report_, allele_status, record, pos);
      continue;
    }
    const Unsupported genotype_status = genotype_reason(
        selected_sample_fields, gt, format_count, line_number);
    if (genotype_status != Unsupported::none) {
      unsupported_record(skip_unsupported, report_, genotype_status, record, pos);
      continue;
    }

    if (chromosomes_.empty() || chromosomes_.back().chromosome != record[0]) {
      chromosomes_.push_back({record[0],
                              static_cast<std::uint64_t>(variants_.size()), 0u});
    }
    variants_.push_back({record[0], id, pos, record[3], record[4], generated});
    record_lines_.push_back(line_number);
    ++chromosomes_.back().variant_count;
    ++report_.retained_variants;
  }
  input->finish();
  report_.maximum_parsing_buffer_bytes = maximum_line_bytes + input->fixed_buffer_bytes();
  if (!header) invalid("missing #CHROM header");
  if (!saw_selected_chromosome) {
    invalid("selected chromosome '" + selected_chromosome + "' is absent from the VCF");
  }
  if (variants_.empty()) invalid("no supported variants were retained");
}

PhasedVcfReader::~PhasedVcfReader() = default;

const std::vector<std::string>& PhasedVcfReader::samples() const noexcept {
  return samples_;
}

const std::vector<VcfVariant>& PhasedVcfReader::variants() const noexcept {
  return variants_;
}

const std::vector<VcfChromosomeBlock>& PhasedVcfReader::chromosomes() const noexcept {
  return chromosomes_;
}

const VcfImportReport& PhasedVcfReader::report() const noexcept {
  return report_;
}

void PhasedVcfReader::start_chromosome(std::uint64_t chromosome_index) {
  if (chromosome_index >= chromosomes_.size()) invalid("chromosome index is out of range");
  const auto& block = chromosomes_[static_cast<std::size_t>(chromosome_index)];
  cursor_ = block.first_variant;
  cursor_end_ = cursor_ + block.variant_count;
  stream_ = open_input(path_, report_.input_type);
  stream_line_ = 0u;
}

bool PhasedVcfReader::next(std::uint64_t& variant_index,
                           std::vector<std::uint8_t>& h1,
                           std::vector<std::uint8_t>& h2) {
  if (cursor_ >= cursor_end_) return false;
  if (!stream_) invalid("chromosome streaming has not been started");
  const std::uint64_t target_line = record_lines_[static_cast<std::size_t>(cursor_)];
  std::string line;
  while (stream_line_ < target_line) {
    if (!stream_->next(line)) io_failure("input ended before a retained VCF record");
    ++stream_line_;
  }
  const auto record = record_fields(
      line, stream_line_, static_cast<std::size_t>(report_.vcf_sample_count));
  const auto& expected = variants_[static_cast<std::size_t>(cursor_)];
  const std::uint64_t pos = position(record[1], stream_line_);
  std::string id = record[2] == "." ?
      record[0] + ":" + record[1] + ":" + record[3] + ":" + record[4] : record[2];
  if (record[0] != expected.chromosome || id != expected.variant_id ||
      pos != expected.base_pair_position || record[3] != expected.reference_allele ||
      record[4] != expected.alternate_allele) {
    io_failure("VCF changed after its validation scan");
  }
  std::size_t format_count = 0;
  const std::size_t gt = gt_index(record[8], stream_line_, format_count);
  const auto selected_sample_fields =
      selected_fields(line, record, selected_columns_, stream_line_);
  h1.assign(samples_.size(), 0u);
  h2.assign(samples_.size(), 0u);
  for (std::size_t sample = 0; sample < selected_columns_.size(); ++sample) {
    const Unsupported status = classify_gt(
        selected_gt(selected_sample_fields[sample], gt, format_count, stream_line_),
        &h1[sample], &h2[sample], stream_line_);
    if (status != Unsupported::none) io_failure("VCF changed after its validation scan");
  }
  variant_index = cursor_++;
  return true;
}

}  // namespace gsim::native::metadata
