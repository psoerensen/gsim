#ifndef GSIM_STANDALONE_METADATA_HPP
#define GSIM_STANDALONE_METADATA_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace gsim::native::metadata {

struct MarkerMetadata {
  std::string chromosome;
  std::string marker_id;
  double genetic_distance_cm;
  std::uint64_t base_pair_position;
  std::string allele1;
  std::string allele2;
};

struct SimulationSampleMetadata {
  std::string family_id;
  std::string individual_id;
  std::string paternal_id;
  std::string maternal_id;
  std::uint32_t sex;
};

struct MetadataWriteResult {
  std::uint64_t record_count;
  std::uint64_t bytes_written;
  std::uint64_t maximum_record_bytes;
};

class ValidatedVariantMetadata final {
 public:
  explicit ValidatedVariantMetadata(std::vector<MarkerMetadata> records);
  [[nodiscard]] static ValidatedVariantMetadata read_bim(
      const std::string& path);

  [[nodiscard]] std::uint64_t size() const noexcept;
  [[nodiscard]] const std::vector<MarkerMetadata>& records() const noexcept;
  [[nodiscard]] MetadataWriteResult write_bim(const std::string& path,
                                               bool overwrite) const;

 private:
  std::vector<MarkerMetadata> records_;
};

class ValidatedSampleMetadata final {
 public:
  explicit ValidatedSampleMetadata(
      std::vector<SimulationSampleMetadata> records);
  [[nodiscard]] static ValidatedSampleMetadata read_fam(
      const std::string& path);

  [[nodiscard]] std::uint64_t size() const noexcept;
  [[nodiscard]] const std::vector<SimulationSampleMetadata>& records()
      const noexcept;
  [[nodiscard]] MetadataWriteResult write_fam(const std::string& path,
                                               bool overwrite) const;

 private:
  std::vector<SimulationSampleMetadata> records_;
};

}  // namespace gsim::native::metadata

#endif
