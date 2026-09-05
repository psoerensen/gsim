#include "standalone_vcf.h"

#include "standalone_metadata_error.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <sstream>
#include <unordered_set>

namespace gsim::native::metadata {
namespace {

namespace fs = std::filesystem;

[[noreturn]] void invalid(const std::string& message) {
  throw Error(StatusCode::invalid_argument, "VCF: " + message);
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

void allele(const std::string& value, const char* name, std::uint64_t line) {
  if (value.size() != 1u ||
      (value[0] != 'A' && value[0] != 'C' && value[0] != 'G' &&
       value[0] != 'T')) {
    invalid(std::string(name) +
            " must be one uppercase A/C/G/T base at line " +
            std::to_string(line));
  }
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
    if (values[i] == "GT") {
      if (found != values.size()) {
        invalid("FORMAT contains duplicate GT at line " +
                std::to_string(line));
      }
      found = i;
    }
  }
  if (found == values.size()) {
    invalid("FORMAT does not contain GT at line " + std::to_string(line));
  }
  return found;
}

void parse_gt(const std::string& value, std::uint8_t& left,
              std::uint8_t& right, std::uint64_t line) {
  if (value.size() != 3u || value[1] != '|' ||
      (value[0] != '0' && value[0] != '1') ||
      (value[2] != '0' && value[2] != '1')) {
    invalid("GT must be phased diploid 0|0, 0|1, 1|0, or 1|1 at line " +
            std::to_string(line));
  }
  left = static_cast<std::uint8_t>(value[0] - '0');
  right = static_cast<std::uint8_t>(value[2] - '0');
}

void parse_samples(const std::vector<std::string>& record,
                   std::size_t sample_count, std::uint64_t line,
                   std::vector<std::uint8_t>* h1,
                   std::vector<std::uint8_t>* h2) {
  if (record.size() != 9u + sample_count) {
    invalid("sample-field count does not match #CHROM header at line " +
            std::to_string(line));
  }
  std::size_t format_count = 0;
  const std::size_t gt = gt_index(record[8], line, format_count);
  if (format_count == 0u) invalid("empty FORMAT at line " + std::to_string(line));
  if (h1 != nullptr) {
    h1->assign(sample_count, 0u);
    h2->assign(sample_count, 0u);
  }
  for (std::size_t sample = 0; sample < sample_count; ++sample) {
    const auto values = fields(record[9u + sample], ':');
    if (values.size() != format_count) {
      invalid("sample FORMAT arity mismatch at line " +
              std::to_string(line));
    }
    std::uint8_t left = 0;
    std::uint8_t right = 0;
    parse_gt(values[gt], left, right, line);
    if (h1 != nullptr) {
      (*h1)[sample] = left;
      (*h2)[sample] = right;
    }
  }
}

std::vector<std::string> record_fields(std::string line, std::uint64_t number) {
  if (!line.empty() && line.back() == '\r') line.pop_back();
  const auto result = fields(line, '\t');
  if (result.size() < 10u) {
    invalid("record has fewer than ten tab-delimited fields at line " +
            std::to_string(number));
  }
  return result;
}

std::uint64_t offset(std::streampos value) {
  if (value < std::streampos(0)) {
    throw Error(StatusCode::internal_error, "VCF stream offset is invalid");
  }
  const auto converted = static_cast<unsigned long long>(value);
  if (converted > std::numeric_limits<std::uint64_t>::max()) {
    throw Error(StatusCode::invalid_extent, "VCF stream offset exceeds uint64_t");
  }
  return static_cast<std::uint64_t>(converted);
}

}  // namespace

PhasedVcfReader::PhasedVcfReader(const std::string& path) : path_(path) {
  if (path.empty()) invalid("path must not be empty");
  std::string lower = path;
  std::transform(lower.begin(), lower.end(), lower.begin(),
                 [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
  if (lower.size() >= 7u && lower.substr(lower.size() - 7u) == ".vcf.gz") {
    invalid(".vcf.gz is unsupported; provide ordinary uncompressed .vcf");
  }
  stream_.open(fs::u8path(path_), std::ios::binary);
  if (!stream_.is_open()) {
    throw Error(StatusCode::internal_error, "could not open VCF: " + path_);
  }
  std::unordered_set<std::string> sample_ids;
  std::unordered_set<std::string> variant_ids;
  std::unordered_set<std::string> closed_chromosomes;
  std::string current_chromosome;
  std::uint64_t previous_position = 0;
  std::uint64_t line_number = 0;
  bool header = false;
  bool fileformat = false;
  std::string line;
  while (true) {
    const std::streampos start = stream_.tellg();
    if (!std::getline(stream_, line)) break;
    ++line_number;
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
      samples_.assign(heading.begin() + 9, heading.end());
      if (!fileformat) invalid("missing ##fileformat=VCFv declaration");
      for (const auto& sample : samples_) {
        if (sample.empty() || sample == "." || !sample_ids.insert(sample).second) {
          invalid("sample IDs must be unique, nonempty, and not '.'");
        }
      }
      header = true;
      continue;
    }
    if (line.empty() || line[0] == '#') {
      invalid("blank or header line after #CHROM at line " +
              std::to_string(line_number));
    }
    const auto record = record_fields(line, line_number);
    if (record[0].empty()) invalid("empty chromosome at line " + std::to_string(line_number));
    const std::uint64_t pos = position(record[1], line_number);
    allele(record[3], "REF", line_number);
    allele(record[4], "ALT", line_number);
    if (record[3] == record[4]) invalid("REF and ALT are identical at line " + std::to_string(line_number));
    parse_samples(record, samples_.size(), line_number, nullptr, nullptr);

    if (record[0] != current_chromosome) {
      if (!current_chromosome.empty()) closed_chromosomes.insert(current_chromosome);
      if (closed_chromosomes.count(record[0]) != 0u) {
        invalid("chromosome '" + record[0] + "' occurs in disjoint blocks");
      }
      current_chromosome = record[0];
      previous_position = 0;
      chromosomes_.push_back({record[0],
                              static_cast<std::uint64_t>(variants_.size()), 0u});
    }
    if (pos < previous_position) {
      invalid("POS decreases within chromosome '" + current_chromosome + "'");
    }
    previous_position = pos;
    std::string id = record[2];
    if (id.empty()) invalid("empty variant ID at line " + std::to_string(line_number));
    const bool generated = id == ".";
    if (generated) {
      id = record[0] + ":" + record[1] + ":" + record[3] + ":" + record[4];
    }
    if (!variant_ids.insert(id).second) invalid("duplicate final variant ID '" + id + "'");
    variants_.push_back({record[0], id, pos, record[3], record[4], generated});
    offsets_.push_back(offset(start));
    ++chromosomes_.back().variant_count;
  }
  if (!header) invalid("missing #CHROM header");
  if (variants_.empty()) invalid("contains no variant records");
  stream_.clear();
}

const std::vector<std::string>& PhasedVcfReader::samples() const noexcept {
  return samples_;
}

const std::vector<VcfVariant>& PhasedVcfReader::variants() const noexcept {
  return variants_;
}

const std::vector<VcfChromosomeBlock>& PhasedVcfReader::chromosomes() const noexcept {
  return chromosomes_;
}

void PhasedVcfReader::start_chromosome(std::uint64_t chromosome_index) {
  if (chromosome_index >= chromosomes_.size()) invalid("chromosome index is out of range");
  const auto& block = chromosomes_[static_cast<std::size_t>(chromosome_index)];
  cursor_ = block.first_variant;
  cursor_end_ = cursor_ + block.variant_count;
}

bool PhasedVcfReader::next(std::uint64_t& variant_index,
                           std::vector<std::uint8_t>& h1,
                           std::vector<std::uint8_t>& h2) {
  if (cursor_ >= cursor_end_) return false;
  stream_.clear();
  const std::uint64_t stored = offsets_[static_cast<std::size_t>(cursor_)];
  if (stored > static_cast<std::uint64_t>(
                   std::numeric_limits<std::streamoff>::max())) {
    throw Error(StatusCode::invalid_extent, "VCF stream offset exceeds streamoff");
  }
  stream_.seekg(static_cast<std::streamoff>(stored));
  std::string line;
  if (!std::getline(stream_, line)) {
    throw Error(StatusCode::internal_error, "could not reread VCF record");
  }
  const auto record = record_fields(line, cursor_ + 2u);
  parse_samples(record, samples_.size(), cursor_ + 2u, &h1, &h2);
  variant_index = cursor_++;
  return true;
}

}  // namespace gsim::native::metadata
