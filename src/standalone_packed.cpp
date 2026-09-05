#include "standalone_packed.h"

#include "standalone_error.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace gsim::native {
namespace {

std::uint64_t checked_extent(std::uint64_t rows, std::uint64_t columns,
                             std::uint64_t leading_dimension,
                             std::uint64_t buffer_length,
                             const char* name) {
    if (rows == 0u || columns == 0u) {
        invalid_argument(std::string(name) + " dimensions must be positive");
    }
    if (leading_dimension < rows) {
        invalid_argument(std::string(name) +
                         " leading dimension is smaller than row count");
    }
    const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
    if (columns - 1u > (maximum - rows) / leading_dimension) {
        invalid_argument(std::string(name) + " extent overflows uint64_t");
    }
    const std::uint64_t extent = (columns - 1u) * leading_dimension + rows;
    if (buffer_length < extent) {
        invalid_argument(std::string(name) + " buffer is too short");
    }
    if (extent > static_cast<std::uint64_t>(
                     std::numeric_limits<std::size_t>::max())) {
        invalid_argument(std::string(name) +
                         " extent exceeds platform indexing limits");
    }
    return extent;
}

void require_compatible(const PhasedHaplotypeMatrix& first,
                        const PhasedHaplotypeMatrix& second,
                        const char* operation) {
    if (first.individual_count() != second.individual_count() ||
        first.marker_count() != second.marker_count()) {
        invalid_argument(std::string(operation) +
                         " requires identical matrix dimensions");
    }
}

} // namespace

PhasedHaplotypeMatrix::PhasedHaplotypeMatrix(std::uint64_t individuals,
                                             std::uint64_t markers)
    : individuals_(individuals), markers_(markers), words_per_marker_(0u) {
    if (individuals == 0u || markers == 0u) {
        invalid_argument("phased haplotype dimensions must be positive");
    }
    words_per_marker_ = individuals / 64u +
                        static_cast<std::uint64_t>(individuals % 64u != 0u);
    const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
    if (markers > maximum / words_per_marker_) {
        invalid_argument("phased haplotype word count overflows uint64_t");
    }
    const std::uint64_t total_words = markers * words_per_marker_;
    if (total_words > static_cast<std::uint64_t>(
                          std::numeric_limits<std::size_t>::max() /
                          sizeof(std::uint64_t))) {
        invalid_argument("phased haplotype allocation exceeds platform limits");
    }
    words_.assign(static_cast<std::size_t>(total_words), 0u);
}

PhasedHaplotypeMatrix PhasedHaplotypeMatrix::from_values(
    std::uint64_t individuals, std::uint64_t markers,
    const std::uint8_t* values, std::uint64_t values_length,
    std::uint64_t leading_dimension) {
    if (values == nullptr) {
        invalid_argument("phased haplotype values must not be null");
    }
    (void)checked_extent(individuals, markers, leading_dimension,
                         values_length, "phased haplotype input");
    PhasedHaplotypeMatrix result(individuals, markers);
    for (std::uint64_t marker = 0; marker < markers; ++marker) {
        for (std::uint64_t individual = 0; individual < individuals;
             ++individual) {
            const std::uint8_t value = values[static_cast<std::size_t>(
                marker * leading_dimension + individual)];
            if (value > 1u) {
                invalid_argument("phased haplotype input contains a nonbinary allele");
            }
            result.set_allele(individual, marker, value);
        }
    }
    return result;
}

std::uint64_t PhasedHaplotypeMatrix::storage_bytes() const noexcept {
    return static_cast<std::uint64_t>(words_.size()) *
           static_cast<std::uint64_t>(sizeof(std::uint64_t));
}

bool PhasedHaplotypeMatrix::has_canonical_padding() const noexcept {
    const unsigned int used = static_cast<unsigned int>(individuals_ % 64u);
    if (used == 0u) return true;
    const std::uint64_t padding_mask =
        ~((std::uint64_t{1} << used) - std::uint64_t{1});
    for (std::uint64_t marker = 0u; marker < markers_; ++marker) {
        const std::size_t last_word = static_cast<std::size_t>(
            marker * words_per_marker_ + words_per_marker_ - 1u);
        if ((words_[last_word] & padding_mask) != 0u) return false;
    }
    return true;
}

std::size_t PhasedHaplotypeMatrix::offset(std::uint64_t individual,
                                         std::uint64_t marker) const noexcept {
    return static_cast<std::size_t>(
        marker * words_per_marker_ + individual / 64u);
}

void PhasedHaplotypeMatrix::validate_cell(std::uint64_t individual,
                                          std::uint64_t marker) const {
    if (individual >= individuals_ || marker >= markers_) {
        throw Error(GBITS_STATUS_OUT_OF_RANGE,
                    "phased haplotype cell index is out of range");
    }
}

void PhasedHaplotypeMatrix::validate_interval(
    std::uint64_t first_marker, std::uint64_t last_marker) const {
    if (first_marker > last_marker || last_marker >= markers_) {
        throw Error(GBITS_STATUS_OUT_OF_RANGE,
                    "phased haplotype interval is invalid");
    }
}

std::uint64_t PhasedHaplotypeMatrix::word(
    std::uint64_t marker, std::uint64_t word_index) const {
    if (marker >= markers_ || word_index >= words_per_marker_) {
        throw Error(GBITS_STATUS_OUT_OF_RANGE,
                    "phased haplotype word index is out of range");
    }
    return words_[static_cast<std::size_t>(marker * words_per_marker_ +
                                           word_index)];
}

void PhasedHaplotypeMatrix::set_word(
    std::uint64_t marker, std::uint64_t word_index, std::uint64_t value) {
    if (marker >= markers_ || word_index >= words_per_marker_) {
        throw Error(GBITS_STATUS_OUT_OF_RANGE,
                    "phased haplotype word index is out of range");
    }
    if (word_index + 1u == words_per_marker_ && individuals_ % 64u != 0u) {
        const unsigned int used = static_cast<unsigned int>(individuals_ % 64u);
        const std::uint64_t padding_mask =
            ~((std::uint64_t{1} << used) - std::uint64_t{1});
        if ((value & padding_mask) != 0u) {
            invalid_argument("phased haplotype word has noncanonical padding");
        }
    }
    words_[static_cast<std::size_t>(marker * words_per_marker_ +
                                    word_index)] = value;
}

std::uint8_t PhasedHaplotypeMatrix::allele(
    std::uint64_t individual, std::uint64_t marker) const {
    validate_cell(individual, marker);
    const unsigned int bit = static_cast<unsigned int>(individual % 64u);
    return static_cast<std::uint8_t>((words_[offset(individual, marker)] >> bit) &
                                     std::uint64_t{1});
}

void PhasedHaplotypeMatrix::set_allele(
    std::uint64_t individual, std::uint64_t marker, std::uint8_t value) {
    validate_cell(individual, marker);
    if (value > 1u) {
        invalid_argument("phased haplotype allele must be zero or one");
    }
    const unsigned int bit = static_cast<unsigned int>(individual % 64u);
    const std::uint64_t mask = std::uint64_t{1} << bit;
    std::uint64_t& target = words_[offset(individual, marker)];
    target = value == 0u ? target & ~mask : target | mask;
}

void PhasedHaplotypeMatrix::unpack(
    std::uint8_t* output, std::uint64_t output_length,
    std::uint64_t leading_dimension) const {
    if (output == nullptr) {
        invalid_argument("phased haplotype output must not be null");
    }
    (void)checked_extent(individuals_, markers_, leading_dimension,
                         output_length, "phased haplotype output");
    for (std::uint64_t marker = 0; marker < markers_; ++marker) {
        for (std::uint64_t individual = 0; individual < individuals_;
             ++individual) {
            output[static_cast<std::size_t>(marker * leading_dimension +
                                            individual)] =
                allele(individual, marker);
        }
    }
}

void PhasedHaplotypeMatrix::copy_interval(
    std::uint64_t destination_individual,
    const PhasedHaplotypeMatrix& source,
    std::uint64_t source_individual,
    std::uint64_t first_marker,
    std::uint64_t last_marker) {
    if (markers_ != source.markers_) {
        invalid_argument("interval copy requires equal marker counts");
    }
    validate_cell(destination_individual, 0u);
    source.validate_cell(source_individual, 0u);
    validate_interval(first_marker, last_marker);
    for (std::uint64_t marker = first_marker;; ++marker) {
        set_allele(destination_individual, marker,
                   source.allele(source_individual, marker));
        if (marker == last_marker) break;
    }
}

void PhasedHaplotypeMatrix::copy_filtered_segment(
    std::uint64_t destination_individual,
    const PhasedHaplotypeMatrix& source,
    std::uint64_t source_individual,
    std::uint64_t first_marker,
    std::uint64_t last_marker,
    double coalescent_age,
    const double* mutation_age,
    std::uint64_t mutation_age_count) {
    if (mutation_age == nullptr || mutation_age_count != markers_) {
        invalid_argument("mutation ages must exactly match marker count");
    }
    if (!std::isfinite(coalescent_age) || coalescent_age < 0.0) {
        invalid_argument("coalescent age must be finite and nonnegative");
    }
    if (markers_ != source.markers_) {
        invalid_argument("filtered copy requires equal marker counts");
    }
    validate_cell(destination_individual, 0u);
    source.validate_cell(source_individual, 0u);
    validate_interval(first_marker, last_marker);
    for (std::uint64_t marker = first_marker;; ++marker) {
        if (!std::isfinite(mutation_age[static_cast<std::size_t>(marker)]) ||
            mutation_age[static_cast<std::size_t>(marker)] < 0.0) {
            invalid_argument("mutation ages must be finite and nonnegative");
        }
        if (marker == last_marker) break;
    }
    for (std::uint64_t marker = first_marker;; ++marker) {
        const std::uint8_t retained =
            source.allele(source_individual, marker) == 1u &&
                    coalescent_age < mutation_age[static_cast<std::size_t>(marker)]
                ? std::uint8_t{1}
                : std::uint8_t{0};
        set_allele(destination_individual, marker, retained);
        if (marker == last_marker) break;
    }
}

void PhasedHaplotypeMatrix::make_gamete(
    std::uint64_t destination_individual,
    const PhasedHaplotypeMatrix& parent_h1,
    const PhasedHaplotypeMatrix& parent_h2,
    std::uint64_t parent_individual,
    std::uint32_t starting_haplotype,
    const std::uint64_t* crossover_boundaries,
    std::uint64_t crossover_count) {
    require_compatible(parent_h1, parent_h2, "gamete construction");
    if (markers_ != parent_h1.markers_) {
        invalid_argument("gamete construction requires equal marker counts");
    }
    validate_cell(destination_individual, 0u);
    parent_h1.validate_cell(parent_individual, 0u);
    if (starting_haplotype > 1u) {
        invalid_argument("starting haplotype must be zero or one");
    }
    if ((crossover_boundaries == nullptr) != (crossover_count == 0u)) {
        invalid_argument("crossover boundary pointer/count combination is invalid");
    }
    std::uint64_t previous = 0u;
    for (std::uint64_t i = 0; i < crossover_count; ++i) {
        const std::uint64_t boundary =
            crossover_boundaries[static_cast<std::size_t>(i)];
        if (boundary > markers_ || (i > 0u && boundary < previous)) {
            invalid_argument("crossover boundaries must be sorted and in range");
        }
        previous = boundary;
    }
    std::uint32_t source = starting_haplotype;
    std::uint64_t crossover = 0u;
    for (std::uint64_t marker = 0; marker < markers_; ++marker) {
        while (crossover < crossover_count &&
               crossover_boundaries[static_cast<std::size_t>(crossover)] ==
                   marker) {
            source ^= 1u;
            ++crossover;
        }
        set_allele(destination_individual, marker,
                   source == 0u ? parent_h1.allele(parent_individual, marker)
                                : parent_h2.allele(parent_individual, marker));
    }
}

void PhasedHaplotypeMatrix::decode_genotypes(
    const PhasedHaplotypeMatrix& h2,
    std::uint8_t* output,
    std::uint64_t output_length,
    std::uint64_t leading_dimension) const {
    require_compatible(*this, h2, "genotype decoding");
    if (output == nullptr) {
        invalid_argument("genotype output must not be null");
    }
    (void)checked_extent(individuals_, markers_, leading_dimension,
                         output_length, "genotype output");
    for (std::uint64_t marker = 0; marker < markers_; ++marker) {
        for (std::uint64_t individual = 0; individual < individuals_;
             ++individual) {
            output[static_cast<std::size_t>(marker * leading_dimension +
                                            individual)] =
                static_cast<std::uint8_t>(allele(individual, marker) +
                                          h2.allele(individual, marker));
        }
    }
}

} // namespace gsim::native
