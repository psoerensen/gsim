#ifndef GSIM_STANDALONE_VCF_HPP
#define GSIM_STANDALONE_VCF_HPP

#include <cstdint>
#include <fstream>
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

class PhasedVcfReader final {
 public:
  explicit PhasedVcfReader(const std::string& path);
  [[nodiscard]] const std::vector<std::string>& samples() const noexcept;
  [[nodiscard]] const std::vector<VcfVariant>& variants() const noexcept;
  [[nodiscard]] const std::vector<VcfChromosomeBlock>& chromosomes() const noexcept;
  void start_chromosome(std::uint64_t chromosome_index);
  [[nodiscard]] bool next(std::uint64_t& variant_index,
                          std::vector<std::uint8_t>& h1,
                          std::vector<std::uint8_t>& h2);

 private:
  std::string path_;
  std::vector<std::string> samples_;
  std::vector<VcfVariant> variants_;
  std::vector<VcfChromosomeBlock> chromosomes_;
  std::vector<std::uint64_t> offsets_;
  std::ifstream stream_;
  std::uint64_t cursor_{0};
  std::uint64_t cursor_end_{0};
};

}  // namespace gsim::native::metadata

#endif
