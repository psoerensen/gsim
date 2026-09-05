#ifndef GSIM_STANDALONE_PACKED_HPP
#define GSIM_STANDALONE_PACKED_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

namespace gsim::native {

class PhasedHaplotypeMatrix final {
public:
    PhasedHaplotypeMatrix(std::uint64_t individuals,
                          std::uint64_t markers);

    static PhasedHaplotypeMatrix from_values(
        std::uint64_t individuals, std::uint64_t markers,
        const std::uint8_t* values, std::uint64_t values_length,
        std::uint64_t leading_dimension);

    std::uint64_t individual_count() const noexcept { return individuals_; }
    std::uint64_t marker_count() const noexcept { return markers_; }
    std::uint64_t words_per_marker() const noexcept { return words_per_marker_; }
    std::uint64_t storage_bytes() const noexcept;
    bool has_canonical_padding() const noexcept;
    std::uint64_t word(std::uint64_t marker,
                       std::uint64_t word_index) const;
    void set_word(std::uint64_t marker, std::uint64_t word_index,
                  std::uint64_t value);
    std::uint8_t allele(std::uint64_t individual,
                        std::uint64_t marker) const;
    void set_allele(std::uint64_t individual, std::uint64_t marker,
                    std::uint8_t allele);
    void unpack(std::uint8_t* output, std::uint64_t output_length,
                std::uint64_t leading_dimension) const;
    void copy_interval(std::uint64_t destination_individual,
                       const PhasedHaplotypeMatrix& source,
                       std::uint64_t source_individual,
                       std::uint64_t first_marker,
                       std::uint64_t last_marker);
    void copy_filtered_segment(std::uint64_t destination_individual,
                               const PhasedHaplotypeMatrix& source,
                               std::uint64_t source_individual,
                               std::uint64_t first_marker,
                               std::uint64_t last_marker,
                               double coalescent_age,
                               const double* mutation_age,
                               std::uint64_t mutation_age_count);
    void make_gamete(std::uint64_t destination_individual,
                     const PhasedHaplotypeMatrix& parent_h1,
                     const PhasedHaplotypeMatrix& parent_h2,
                     std::uint64_t parent_individual,
                     std::uint32_t starting_haplotype,
                     const std::uint64_t* crossover_boundaries,
                     std::uint64_t crossover_count);
    void decode_genotypes(const PhasedHaplotypeMatrix& h2,
                          std::uint8_t* output,
                          std::uint64_t output_length,
                          std::uint64_t leading_dimension) const;

private:
    std::size_t offset(std::uint64_t individual,
                       std::uint64_t marker) const noexcept;
    void validate_cell(std::uint64_t individual,
                       std::uint64_t marker) const;
    void validate_interval(std::uint64_t first_marker,
                           std::uint64_t last_marker) const;

    std::uint64_t individuals_;
    std::uint64_t markers_;
    std::uint64_t words_per_marker_;
    std::vector<std::uint64_t> words_;
};

} // namespace gsim::native

#endif
