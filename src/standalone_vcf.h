#ifndef GSIM_STANDALONE_VCF_H
#define GSIM_STANDALONE_VCF_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace gsim::native::metadata {

struct VcfVariant {
  std::string chromosome;
  std::string variant_id;
  std::uint64_t base_pair_position;
  std::string reference_allele;
  std::string alternate_allele;
  bool generated_id;
};

struct VcfChromosomeBlock {
  std::string chromosome;
  std::uint64_t first_variant;
  std::uint64_t variant_count;
};

struct VcfImportReport {
  std::string input_type;
  std::uint64_t vcf_sample_count{0};
  std::uint64_t selected_sample_count{0};
  std::uint64_t total_records_scanned{0};
  std::uint64_t retained_variants{0};
  std::uint64_t outside_selected_chromosome{0};
  std::uint64_t outside_selected_region{0};
  std::uint64_t indels{0};
  std::uint64_t multiallelic_records{0};
  std::uint64_t symbolic_or_breakend_alleles{0};
  std::uint64_t other_unsupported_alleles{0};
  std::uint64_t missing_gt{0};
  std::uint64_t unphased_gt{0};
  std::uint64_t non_diploid_gt{0};
  std::uint64_t duplicate_final_ids{0};
  std::uint64_t maximum_parsing_buffer_bytes{0};
};

class VcfLineInput;

class PhasedVcfReader final {
 public:
  PhasedVcfReader(const std::string& path,
                  const std::vector<std::string>& selected_samples,
                  const std::string& selected_chromosome,
                  bool has_region, std::uint64_t region_start,
                  std::uint64_t region_end, bool skip_unsupported);
  ~PhasedVcfReader();

  PhasedVcfReader(const PhasedVcfReader&) = delete;
  PhasedVcfReader& operator=(const PhasedVcfReader&) = delete;

  [[nodiscard]] const std::vector<std::string>& samples() const noexcept;
  [[nodiscard]] const std::vector<VcfVariant>& variants() const noexcept;
  [[nodiscard]] const std::vector<VcfChromosomeBlock>& chromosomes() const noexcept;
  [[nodiscard]] const VcfImportReport& report() const noexcept;
  void start_chromosome(std::uint64_t chromosome_index);
  [[nodiscard]] bool next(std::uint64_t& variant_index,
                          std::vector<std::uint8_t>& h1,
                          std::vector<std::uint8_t>& h2);

 private:
  std::string path_;
  std::vector<std::string> samples_;
  std::vector<std::size_t> selected_columns_;
  std::vector<VcfVariant> variants_;
  std::vector<VcfChromosomeBlock> chromosomes_;
  std::vector<std::uint64_t> record_lines_;
  VcfImportReport report_;
  std::unique_ptr<VcfLineInput> stream_;
  std::uint64_t stream_line_{0};
  std::uint64_t cursor_{0};
  std::uint64_t cursor_end_{0};
};

}  // namespace gsim::native::metadata

#endif
