#include "standalone_metadata.h"

#include "standalone_metadata_error.h"

#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace gsim::native::metadata {
namespace {

namespace fs = std::filesystem;

std::atomic<std::uint64_t> temporary_sequence{0};

std::uint64_t checked_size(std::size_t value, const char* resource) {
  if (value > std::numeric_limits<std::uint64_t>::max()) {
    throw Error(StatusCode::invalid_extent,
                std::string(resource) + " count exceeds uint64_t");
  }
  return static_cast<std::uint64_t>(value);
}

void require_field(const std::string& value, const char* field) {
  if (value.empty()) {
    throw Error(StatusCode::invalid_argument,
                std::string(field) + " must not be empty");
  }
  for (const unsigned char byte : value) {
    if (byte <= 0x20u || byte == 0x7fu) {
      throw Error(StatusCode::invalid_argument,
                  std::string(field) +
                      " must not contain ASCII whitespace or control bytes");
    }
  }
}

void require_allele(const std::string& value, const char* field) {
  require_field(value, field);
  if (value.size() != 1u ||
      (value[0] != 'A' && value[0] != 'C' && value[0] != 'G' &&
       value[0] != 'T')) {
    throw Error(StatusCode::invalid_argument,
                std::string(field) +
                    " must be one uppercase single-base allele A, C, G, or T");
  }
}

std::string format_cm(double value) {
  if (value == 0.0) return "0";
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << std::fixed << std::setprecision(15) << value;
  std::string result = output.str();
  while (!result.empty() && result.back() == '0') result.pop_back();
  if (!result.empty() && result.back() == '.') result.pop_back();
  return result.empty() ? "0" : result;
}

fs::path utf8_path(const std::string& path) {
  if (path.empty()) {
    throw Error(StatusCode::invalid_argument,
                "metadata destination must not be empty");
  }
  return fs::u8path(path);
}

fs::path temporary_path_for(const fs::path& destination) {
  const auto parent = destination.parent_path();
  if (parent.empty() || !fs::exists(parent) || !fs::is_directory(parent)) {
    throw Error(StatusCode::invalid_argument,
                "metadata destination parent directory does not exist: " +
                    destination.u8string());
  }
  for (unsigned int attempt = 0; attempt < 128u; ++attempt) {
    const std::uint64_t token = ++temporary_sequence;
    fs::path candidate = destination;
    candidate += ".gsim.tmp." + std::to_string(token);
    if (!fs::exists(candidate)) return candidate;
  }
  throw Error(StatusCode::internal_error,
              "could not allocate a unique metadata temporary path: " +
                  destination.u8string());
}

template <typename RecordWriter>
MetadataWriteResult write_records(const std::string& path, bool overwrite,
                                  std::uint64_t count,
                                  RecordWriter&& write_record) {
  const fs::path destination = utf8_path(path);
  if (fs::exists(destination) && !overwrite) {
    throw Error(StatusCode::invalid_argument,
                "metadata destination exists and overwrite is disabled: " +
                    destination.u8string());
  }
  const fs::path temporary = temporary_path_for(destination);
  MetadataWriteResult result{count, 0u, 0u};
  try {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output.is_open()) {
      throw Error(StatusCode::internal_error,
                  "could not create metadata temporary file: " +
                      temporary.u8string());
    }
    for (std::uint64_t index = 0; index < count; ++index) {
      const std::string record = write_record(index);
      if (record.size() >
          static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max())) {
        throw Error(StatusCode::invalid_extent,
                    "metadata record exceeds stream-size range");
      }
      output.write(record.data(), static_cast<std::streamsize>(record.size()));
      if (!output) {
        throw Error(StatusCode::internal_error,
                    "failed while writing metadata temporary file: " +
                        temporary.u8string());
      }
      const std::uint64_t bytes = checked_size(record.size(), "metadata record");
      if (result.bytes_written >
          std::numeric_limits<std::uint64_t>::max() - bytes) {
        throw Error(StatusCode::invalid_extent,
                    "metadata byte count exceeds uint64_t");
      }
      result.bytes_written += bytes;
      if (bytes > result.maximum_record_bytes) {
        result.maximum_record_bytes = bytes;
      }
    }
    output.close();
    if (!output) {
      throw Error(StatusCode::internal_error,
                  "failed while closing metadata temporary file: " +
                      temporary.u8string());
    }
    if (overwrite && fs::exists(destination)) {
      std::error_code remove_error;
      fs::remove(destination, remove_error);
      if (remove_error) {
        throw Error(StatusCode::internal_error,
                    "could not replace metadata destination: " +
                        destination.u8string());
      }
    }
    std::error_code rename_error;
    fs::rename(temporary, destination, rename_error);
    if (rename_error) {
      throw Error(StatusCode::internal_error,
                  "could not publish metadata destination: " +
                      destination.u8string());
    }
    return result;
  } catch (...) {
    std::error_code ignored;
    fs::remove(temporary, ignored);
    throw;
  }
}

std::vector<std::string> split_six_tabs(const std::string& line,
                                        const char* format,
                                        std::uint64_t line_number) {
  std::vector<std::string> fields;
  std::size_t start = 0u;
  for (;;) {
    const std::size_t tab = line.find('\t', start);
    fields.push_back(line.substr(start, tab == std::string::npos
                                           ? std::string::npos
                                           : tab - start));
    if (tab == std::string::npos) break;
    start = tab + 1u;
  }
  if (fields.size() != 6u) {
    throw Error(StatusCode::invalid_argument,
                std::string(format) + " line " + std::to_string(line_number) +
                    " must contain exactly six tab-delimited fields");
  }
  return fields;
}

double parse_cm(const std::string& text, std::uint64_t line_number) {
  std::istringstream input(text);
  input.imbue(std::locale::classic());
  double value = 0.0;
  input >> std::noskipws >> value;
  if (!input || !input.eof() || !std::isfinite(value) || value < 0.0) {
    throw Error(StatusCode::invalid_argument,
                "BIM line " + std::to_string(line_number) +
                    " has an invalid genetic position");
  }
  return value;
}

std::uint64_t parse_u64(const std::string& text, const char* field,
                        std::uint64_t line_number, bool allow_zero) {
  try {
    if (text.empty()) {
      throw std::invalid_argument("invalid unsigned integer");
    }
    for (const unsigned char byte : text) {
      if (byte < static_cast<unsigned char>('0') ||
          byte > static_cast<unsigned char>('9')) {
        throw std::invalid_argument("invalid unsigned integer");
      }
    }
    std::size_t parsed = 0u;
    const unsigned long long value = std::stoull(text, &parsed, 10);
    if (parsed != text.size() || (!allow_zero && value == 0u)) {
      throw std::invalid_argument("invalid unsigned integer");
    }
    return static_cast<std::uint64_t>(value);
  } catch (const std::exception&) {
    throw Error(StatusCode::invalid_argument,
                std::string(field) + " on line " +
                    std::to_string(line_number) + " is invalid");
  }
}

template <typename Record, typename Parser>
std::vector<Record> read_records(const std::string& path, const char* format,
                                 Parser&& parser) {
  if (path.empty()) {
    throw Error(StatusCode::invalid_argument,
                std::string(format) + " source must not be empty");
  }
  std::ifstream input(fs::u8path(path), std::ios::binary);
  if (!input.is_open()) {
    throw Error(StatusCode::invalid_argument,
                "could not open " + std::string(format) + " file: " + path);
  }
  std::vector<Record> records;
  std::string line;
  std::uint64_t line_number = 0u;
  while (std::getline(input, line)) {
    ++line_number;
    records.push_back(parser(split_six_tabs(line, format, line_number),
                             line_number));
  }
  if (!input.eof()) {
    throw Error(StatusCode::invalid_argument,
                "failed while reading " + std::string(format) + " file: " +
                    path);
  }
  return records;
}

}  // namespace

ValidatedVariantMetadata::ValidatedVariantMetadata(
    std::vector<MarkerMetadata> records)
    : records_(std::move(records)) {
  if (records_.empty()) {
    throw Error(StatusCode::invalid_extent,
                "variant metadata must contain at least one record");
  }
  std::unordered_set<std::string> ids;
  std::unordered_set<std::string> completed_chromosomes;
  std::string chromosome;
  double previous_cm = 0.0;
  std::uint64_t previous_bp = 0u;
  for (std::size_t index = 0; index < records_.size(); ++index) {
    const auto& record = records_[index];
    require_field(record.chromosome, "chromosome label");
    require_field(record.marker_id, "variant ID");
    if (!ids.insert(record.marker_id).second) {
      throw Error(StatusCode::invalid_argument,
                  "variant IDs must be globally unique");
    }
    if (!std::isfinite(record.genetic_distance_cm) ||
        record.genetic_distance_cm < 0.0) {
      throw Error(StatusCode::invalid_argument,
                  "genetic positions must be finite and nonnegative");
    }
    if (record.base_pair_position == 0u) {
      throw Error(StatusCode::invalid_argument,
                  "base-pair positions must be positive integers");
    }
    require_allele(record.allele1, "alternate allele");
    require_allele(record.allele2, "reference allele");
    if (record.allele1 == record.allele2) {
      throw Error(StatusCode::invalid_argument,
                  "alternate and reference alleles must differ");
    }
    if (record.chromosome != chromosome) {
      if (!chromosome.empty()) completed_chromosomes.insert(chromosome);
      if (completed_chromosomes.count(record.chromosome) != 0u) {
        throw Error(StatusCode::invalid_argument,
                    "each chromosome label must occupy one contiguous block");
      }
      chromosome = record.chromosome;
      previous_cm = record.genetic_distance_cm;
      previous_bp = record.base_pair_position;
    } else {
      if (record.genetic_distance_cm < previous_cm ||
          record.base_pair_position < previous_bp) {
        throw Error(StatusCode::invalid_argument,
                    "genetic and base-pair positions must be nondecreasing "
                    "within chromosome");
      }
      previous_cm = record.genetic_distance_cm;
      previous_bp = record.base_pair_position;
    }
  }
}

ValidatedVariantMetadata ValidatedVariantMetadata::read_bim(
    const std::string& path) {
  return ValidatedVariantMetadata(read_records<MarkerMetadata>(
      path, "BIM", [](const std::vector<std::string>& fields,
                      std::uint64_t line) {
        return MarkerMetadata{fields[0], fields[1], parse_cm(fields[2], line),
                              parse_u64(fields[3], "base-pair position", line,
                                        false),
                              fields[4], fields[5]};
      }));
}

std::uint64_t ValidatedVariantMetadata::size() const noexcept {
  return static_cast<std::uint64_t>(records_.size());
}

const std::vector<MarkerMetadata>& ValidatedVariantMetadata::records()
    const noexcept {
  return records_;
}

MetadataWriteResult ValidatedVariantMetadata::write_bim(
    const std::string& path, bool overwrite) const {
  return write_records(path, overwrite, size(), [this](std::uint64_t index) {
    const auto& value = records_[static_cast<std::size_t>(index)];
    return value.chromosome + '\t' + value.marker_id + '\t' +
           format_cm(value.genetic_distance_cm) + '\t' +
           std::to_string(value.base_pair_position) + '\t' + value.allele1 +
           '\t' + value.allele2 + '\n';
  });
}

ValidatedSampleMetadata::ValidatedSampleMetadata(
    std::vector<SimulationSampleMetadata> records)
    : records_(std::move(records)) {
  if (records_.empty()) {
    throw Error(StatusCode::invalid_extent,
                "sample metadata must contain at least one record");
  }
  std::unordered_map<std::string, std::size_t> positions;
  for (std::size_t index = 0; index < records_.size(); ++index) {
    const auto& record = records_[index];
    require_field(record.family_id, "family ID");
    require_field(record.individual_id, "individual ID");
    require_field(record.paternal_id, "paternal ID");
    require_field(record.maternal_id, "maternal ID");
    if (!positions.emplace(record.individual_id, index).second) {
      throw Error(StatusCode::invalid_argument,
                  "individual IDs must be globally unique");
    }
    if (record.sex > 2u) {
      throw Error(StatusCode::invalid_argument,
                  "sex must be 0 unknown, 1 male, or 2 female");
    }
  }
  std::unordered_set<std::string> paternal_roles;
  std::unordered_set<std::string> maternal_roles;
  for (std::size_t index = 0; index < records_.size(); ++index) {
    const auto& record = records_[index];
    const bool paternal_missing = record.paternal_id == "0";
    const bool maternal_missing = record.maternal_id == "0";
    if (paternal_missing != maternal_missing) {
      throw Error(StatusCode::invalid_argument,
                  "one-known-parent sample records are not supported");
    }
    if (paternal_missing) continue;
    if (record.paternal_id == record.individual_id ||
        record.maternal_id == record.individual_id) {
      throw Error(StatusCode::invalid_argument,
                  "self-parenting is not allowed");
    }
    if (record.paternal_id == record.maternal_id) {
      throw Error(StatusCode::invalid_argument,
                  "paternal and maternal IDs must differ");
    }
    const auto sire = positions.find(record.paternal_id);
    const auto dam = positions.find(record.maternal_id);
    if (sire == positions.end() || dam == positions.end()) {
      throw Error(StatusCode::invalid_argument,
                  "known parents must be present in the written sample set");
    }
    if (sire->second >= index || dam->second >= index) {
      throw Error(StatusCode::invalid_argument,
                  "known parents must precede offspring in FAM/BED order");
    }
    if (records_[sire->second].sex == 2u) {
      throw Error(StatusCode::invalid_argument,
                  "a known female cannot be recorded as a paternal parent");
    }
    if (records_[dam->second].sex == 1u) {
      throw Error(StatusCode::invalid_argument,
                  "a known male cannot be recorded as a maternal parent");
    }
    paternal_roles.insert(record.paternal_id);
    maternal_roles.insert(record.maternal_id);
  }
  for (const auto& id : paternal_roles) {
    if (maternal_roles.count(id) != 0u) {
      throw Error(StatusCode::invalid_argument,
                  "one individual cannot occupy both parental roles");
    }
  }
}

ValidatedSampleMetadata ValidatedSampleMetadata::read_fam(
    const std::string& path) {
  return ValidatedSampleMetadata(read_records<SimulationSampleMetadata>(
      path, "FAM", [](const std::vector<std::string>& fields,
                      std::uint64_t line) {
        if (fields[5] != "-9") {
          throw Error(StatusCode::invalid_argument,
                      "FAM line " + std::to_string(line) +
                          " must use missing phenotype -9");
        }
        const std::uint64_t sex = parse_u64(fields[4], "sex", line, true);
        if (sex > 2u) {
          throw Error(StatusCode::invalid_argument,
                      "FAM sex must be 0, 1, or 2");
        }
        return SimulationSampleMetadata{
            fields[0], fields[1], fields[2], fields[3],
            static_cast<std::uint32_t>(sex)};
      }));
}

std::uint64_t ValidatedSampleMetadata::size() const noexcept {
  return static_cast<std::uint64_t>(records_.size());
}

const std::vector<SimulationSampleMetadata>&
ValidatedSampleMetadata::records() const noexcept {
  return records_;
}

MetadataWriteResult ValidatedSampleMetadata::write_fam(
    const std::string& path, bool overwrite) const {
  return write_records(path, overwrite, size(), [this](std::uint64_t index) {
    const auto& value = records_[static_cast<std::size_t>(index)];
    return value.family_id + '\t' + value.individual_id + '\t' +
           value.paternal_id + '\t' + value.maternal_id + '\t' +
           std::to_string(value.sex) + "\t-9\n";
  });
}

}  // namespace gsim::native::metadata
