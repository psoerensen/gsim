#ifndef GSIM_PACKED_CHROMOSOME_H
#define GSIM_PACKED_CHROMOSOME_H

#include <cstddef>
#include <cstdint>
#include <utility>
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
    std::pair<std::uint64_t, std::uint64_t> copy_filtered_segment_counts(
        std::uint64_t destination_individual,
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

void materialize_founders(
    PhasedHaplotypeMatrix& destination_h1,
    PhasedHaplotypeMatrix& destination_h2,
    const PhasedHaplotypeMatrix& reference_h1,
    const PhasedHaplotypeMatrix& reference_h2,
    const std::uint64_t* destination, const std::uint32_t* phase,
    const std::uint64_t* donor, const std::uint64_t* first,
    const std::uint64_t* last, const double* age, std::uint64_t event_count,
    const double* mutation, std::uint64_t mutation_count,
    std::uint32_t requested_threads, std::uint64_t* copied,
    std::uint64_t* retained);

} // namespace gsim::native

#endif
