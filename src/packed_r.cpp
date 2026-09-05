#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "bed_storage.h"
#include "hap_storage.h"
#include "packed_chromosome.h"

#include "native_r.h"

namespace {

struct Packed {
    explicit Packed(gsim::native::PhasedHaplotypeMatrix matrix)
        : value(std::move(matrix)) {}
    gsim::native::PhasedHaplotypeMatrix value;
};

struct BedSink {
    BedSink(std::string path, std::uint64_t individuals, bool overwrite,
            std::uint64_t capacity)
        : value(std::move(path), individuals, overwrite, capacity) {}
    gsim::native::PhasedBedSink value;
};

struct HapSink {
    HapSink(std::string path, std::uint64_t individuals, bool overwrite)
        : value(std::move(path), individuals, overwrite) {}
    gsim::native::PhasedHapSink value;
};
struct HapReader {
    explicit HapReader(std::string path) : value(std::move(path)) {}
    gsim::native::PhasedHapReader value;
};

struct BedReader {
    BedReader(const std::string& path, std::uint64_t individuals,
              std::uint64_t markers)
        : input(path, std::ios::binary), n(individuals), m(markers),
          bytes_per_variant((individuals + 3u) / 4u) {
        if (!input || n == 0u || m == 0u) {
            throw std::runtime_error("cannot open bounded BED input");
        }
        std::uint8_t header[3]{};
        input.read(reinterpret_cast<char*>(header), 3);
        if (!input || header[0] != 0x6cu || header[1] != 0x1bu ||
            header[2] != 0x01u) {
            throw std::runtime_error("bounded BED input has invalid header");
        }
    }
    std::ifstream input;
    std::uint64_t n;
    std::uint64_t m;
    std::uint64_t bytes_per_variant;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

Packed* require_packed(SEXP pointer) {
    return gsim::native::r::require<Packed>(pointer, "packed haplotypes");
}

BedSink* require_bed_sink(SEXP pointer) {
    return gsim::native::r::require<BedSink>(pointer, "BED sink");
}

HapSink* require_hap_sink(SEXP pointer) {
    return gsim::native::r::require<HapSink>(pointer, "HAP sink");
}

HapReader* require_hap_reader(SEXP pointer) {
    return gsim::native::r::require<HapReader>(pointer, "HAP reader");
}

int scalar_int(SEXP value, const char* name, int lower = 0) {
    if (TYPEOF(value) != INTSXP || XLENGTH(value) != 1 ||
        INTEGER(value)[0] == NA_INTEGER || INTEGER(value)[0] < lower) {
        fail(std::string(name) + " must be an integer scalar");
    }
    return INTEGER(value)[0];
}

bool scalar_bool(SEXP value, const char* name) {
    if (TYPEOF(value) != LGLSXP || XLENGTH(value) != 1 ||
        LOGICAL(value)[0] == NA_LOGICAL) {
        fail(std::string(name) + " must be one logical value");
    }
    return LOGICAL(value)[0] == TRUE;
}

std::string scalar_utf8(SEXP value, const char* name) {
    if (TYPEOF(value) != STRSXP || XLENGTH(value) != 1 ||
        STRING_ELT(value, 0) == NA_STRING) {
        fail(std::string(name) + " must be one nonmissing string");
    }
    const char* text = Rf_translateCharUTF8(STRING_ELT(value, 0));
    if (text == nullptr || text[0] == '\0') {
        fail(std::string(name) + " must not be empty");
    }
    return text;
}

SEXP make_packed(gsim::native::PhasedHaplotypeMatrix value) {
    return gsim::native::r::make_owned(new Packed(std::move(value)));
}

void set_names(SEXP value, const std::vector<const char*>& names) {
    SEXP nms = PROTECT(Rf_allocVector(STRSXP,
                                      static_cast<R_xlen_t>(names.size())));
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(names.size()); ++i) {
        SET_STRING_ELT(nms, i,
                       Rf_mkChar(names[static_cast<std::size_t>(i)]));
    }
    Rf_setAttrib(value, R_NamesSymbol, nms);
    UNPROTECT(1);
}

double exact_r_number(std::uint64_t value, const char* field) {
    if (value > 9007199254740991ULL) {
        fail(std::string(field) + " exceeds exact R numeric representation");
    }
    return static_cast<double>(value);
}

} // namespace

extern "C" SEXP C_gsim_packed_pack(SEXP values) {
    try {
        if (TYPEOF(values) != RAWSXP || !Rf_isMatrix(values)) {
            fail("packed input must be a raw matrix");
        }
        SEXP dimensions = Rf_getAttrib(values, R_DimSymbol);
        const int individuals = INTEGER(dimensions)[0];
        const int markers = INTEGER(dimensions)[1];
        if (individuals <= 0 || markers <= 0) fail("packed input is empty");
        return make_packed(gsim::native::PhasedHaplotypeMatrix::from_values(
            static_cast<std::uint64_t>(individuals),
            static_cast<std::uint64_t>(markers), RAW(values),
            static_cast<std::uint64_t>(XLENGTH(values)),
            static_cast<std::uint64_t>(individuals)));
    } catch (const std::exception& ex) {
        Rf_error("native packed conversion: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_zero(SEXP individuals_sexp,
                                   SEXP markers_sexp) {
    try {
        const int individuals = scalar_int(individuals_sexp, "individuals", 1);
        const int markers = scalar_int(markers_sexp, "markers", 1);
        return make_packed(gsim::native::PhasedHaplotypeMatrix(
            static_cast<std::uint64_t>(individuals),
            static_cast<std::uint64_t>(markers)));
    } catch (const std::exception& ex) {
        Rf_error("native packed allocation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_set_marker(SEXP h1_pointer, SEXP h2_pointer,
                                          SEXP marker_sexp, SEXP h1_values,
                                          SEXP h2_values) {
    try {
        Packed* h1 = require_packed(h1_pointer);
        Packed* h2 = require_packed(h2_pointer);
        const int marker = scalar_int(marker_sexp, "marker", 1);
        const std::uint64_t individuals = h1->value.individual_count();
        const std::uint64_t markers = h1->value.marker_count();
        if (individuals != h2->value.individual_count() ||
            markers != h2->value.marker_count() ||
            static_cast<std::uint64_t>(marker) > markers) {
            fail("packed phase dimensions or marker index are inconsistent");
        }
        if (TYPEOF(h1_values) != RAWSXP || TYPEOF(h2_values) != RAWSXP ||
            static_cast<std::uint64_t>(XLENGTH(h1_values)) != individuals ||
            static_cast<std::uint64_t>(XLENGTH(h2_values)) != individuals) {
            fail("marker allele buffers must be raw vectors matching individuals");
        }
        const std::uint64_t marker_index = static_cast<std::uint64_t>(marker - 1);
        for (std::uint64_t individual = 0; individual < individuals; ++individual) {
            const std::uint8_t first = RAW(h1_values)[static_cast<R_xlen_t>(individual)];
            const std::uint8_t second = RAW(h2_values)[static_cast<R_xlen_t>(individual)];
            if (first > 1u || second > 1u) fail("marker allele buffers must contain only 0 or 1");
            h1->value.set_allele(individual, marker_index, first);
            h2->value.set_allele(individual, marker_index, second);
        }
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("native packed marker write: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_unpack(SEXP pointer) {
    try {
        Packed* packed = require_packed(pointer);
        const std::uint64_t individuals = packed->value.individual_count();
        const std::uint64_t markers = packed->value.marker_count();
        if (individuals > static_cast<std::uint64_t>(
                                      std::numeric_limits<int>::max()) ||
            markers > static_cast<std::uint64_t>(
                                  std::numeric_limits<int>::max())) {
            fail("packed dimensions exceed R matrix limits");
        }
        SEXP out = PROTECT(Rf_allocMatrix(
            RAWSXP, static_cast<int>(individuals), static_cast<int>(markers)));
        packed->value.unpack(RAW(out),
                             static_cast<std::uint64_t>(XLENGTH(out)),
                             individuals);
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("native bounded unpack: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_info(SEXP pointer) {
    try {
        Packed* packed = require_packed(pointer);
        const std::uint64_t words = packed->value.words_per_marker();
        const std::uint64_t bytes = packed->value.storage_bytes();
        SEXP out = PROTECT(Rf_allocVector(REALSXP, 4));
        REAL(out)[0] = static_cast<double>(packed->value.individual_count());
        REAL(out)[1] = static_cast<double>(packed->value.marker_count());
        REAL(out)[2] = static_cast<double>(words);
        REAL(out)[3] = static_cast<double>(bytes);
        set_names(out, {"individuals", "markers", "words_per_marker",
                        "storage_bytes"});
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("native packed info: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_close(SEXP pointer) {
    try {
        gsim::native::r::release<Packed>(pointer, "packed haplotypes");
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("native packed haplotype close: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_word(SEXP pointer, SEXP marker_sexp,
                                   SEXP word_sexp) {
    try {
        Packed* packed = require_packed(pointer);
        const int marker = scalar_int(marker_sexp, "marker");
        const int word_index = scalar_int(word_sexp, "word_index");
        const std::uint64_t value = packed->value.word(
            static_cast<std::uint64_t>(marker),
            static_cast<std::uint64_t>(word_index));
        SEXP out = PROTECT(Rf_allocVector(RAWSXP, 8));
        for (unsigned int byte = 0; byte < 8u; ++byte) {
            RAW(out)[byte] = static_cast<Rbyte>((value >> (byte * 8u)) & 0xffu);
        }
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("native packed word query: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_copy_interval(
    SEXP destination_pointer, SEXP destination_individual_sexp,
    SEXP source_pointer, SEXP source_individual_sexp,
    SEXP first_marker_sexp, SEXP last_marker_sexp) {
    try {
        Packed* destination = require_packed(destination_pointer);
        Packed* source = require_packed(source_pointer);
        destination->value.copy_interval(
            static_cast<std::uint64_t>(scalar_int(
                destination_individual_sexp, "destination individual")),
            source->value,
            static_cast<std::uint64_t>(scalar_int(
                source_individual_sexp, "source individual")),
            static_cast<std::uint64_t>(scalar_int(first_marker_sexp,
                                                  "first marker")),
            static_cast<std::uint64_t>(scalar_int(last_marker_sexp,
                                                  "last marker")));
        return destination_pointer;
    } catch (const std::exception& ex) {
        Rf_error("native packed interval copy: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_materialize_founders(
    SEXP destination_h1_pointer, SEXP destination_h2_pointer,
    SEXP reference_h1_pointer, SEXP reference_h2_pointer,
    SEXP individuals, SEXP phases, SEXP donors, SEXP starts, SEXP ends,
    SEXP ages, SEXP mutation_age, SEXP individual_offset_sexp,
    SEXP threads_sexp, SEXP return_counts_sexp) {
    try {
        Packed* destination_h1 = require_packed(destination_h1_pointer);
        Packed* destination_h2 = require_packed(destination_h2_pointer);
        Packed* reference_h1 = require_packed(reference_h1_pointer);
        Packed* reference_h2 = require_packed(reference_h2_pointer);
        if ((TYPEOF(individuals) != INTSXP && TYPEOF(individuals) != REALSXP) ||
            TYPEOF(phases) != INTSXP ||
            TYPEOF(donors) != INTSXP || TYPEOF(starts) != INTSXP ||
            TYPEOF(ends) != INTSXP || TYPEOF(ages) != REALSXP ||
            TYPEOF(mutation_age) != REALSXP) {
            fail("founder event columns have invalid storage types");
        }
        const R_xlen_t count = XLENGTH(individuals);
        if (XLENGTH(phases) != count || XLENGTH(donors) != count ||
            XLENGTH(starts) != count || XLENGTH(ends) != count ||
            XLENGTH(ages) != count) {
            fail("founder event columns have inconsistent lengths");
        }
        const int offset = scalar_int(individual_offset_sexp,
                                      "individual_offset");
        const int threads = scalar_int(threads_sexp, "threads", 1);
        const bool return_counts = scalar_bool(return_counts_sexp,
                                               "return_counts");
        std::vector<std::uint64_t> destination(static_cast<std::size_t>(count));
        std::vector<std::uint32_t> phase(static_cast<std::size_t>(count));
        std::vector<std::uint64_t> donor(static_cast<std::size_t>(count));
        std::vector<std::uint64_t> first(static_cast<std::size_t>(count));
        std::vector<std::uint64_t> last(static_cast<std::size_t>(count));
        for (R_xlen_t i = 0; i < count; ++i) {
            const double individual_value = TYPEOF(individuals) == INTSXP
                ? static_cast<double>(INTEGER(individuals)[i])
                : REAL(individuals)[i];
            const int phase_value = INTEGER(phases)[i];
            const int donor_value = INTEGER(donors)[i];
            const int first_value = INTEGER(starts)[i];
            const int last_value = INTEGER(ends)[i];
            if (!R_FINITE(individual_value) ||
                individual_value != std::floor(individual_value) ||
                individual_value <= static_cast<double>(offset) ||
                individual_value > static_cast<double>(std::numeric_limits<int>::max()) ||
                phase_value < 1 || phase_value > 2 ||
                donor_value == NA_INTEGER || donor_value < 1 ||
                first_value == NA_INTEGER || first_value < 1 ||
                last_value == NA_INTEGER || last_value < first_value) {
                fail("founder event plan contains an invalid index");
            }
            destination[static_cast<std::size_t>(i)] =
                static_cast<std::uint64_t>(individual_value) -
                static_cast<std::uint64_t>(offset) - 1u;
            phase[static_cast<std::size_t>(i)] =
                static_cast<std::uint32_t>(phase_value - 1);
            donor[static_cast<std::size_t>(i)] =
                static_cast<std::uint64_t>(donor_value - 1);
            first[static_cast<std::size_t>(i)] =
                static_cast<std::uint64_t>(first_value - 1);
            last[static_cast<std::size_t>(i)] =
                static_cast<std::uint64_t>(last_value - 1);
        }
        std::vector<std::uint64_t> copied;
        std::vector<std::uint64_t> retained;
        if (return_counts) {
            copied.resize(static_cast<std::size_t>(count));
            retained.resize(static_cast<std::size_t>(count));
        }
        gsim::native::materialize_founders(
            destination_h1->value, destination_h2->value,
            reference_h1->value, reference_h2->value,
            destination.empty() ? nullptr : destination.data(),
            phase.empty() ? nullptr : phase.data(),
            donor.empty() ? nullptr : donor.data(),
            first.empty() ? nullptr : first.data(),
            last.empty() ? nullptr : last.data(), REAL(ages),
            static_cast<std::uint64_t>(count), REAL(mutation_age),
            static_cast<std::uint64_t>(XLENGTH(mutation_age)),
            static_cast<std::uint32_t>(threads),
            return_counts ? copied.data() : nullptr,
            return_counts ? retained.data() : nullptr);
        if (!return_counts) return R_NilValue;
        SEXP first_counts = PROTECT(Rf_allocVector(INTSXP, count));
        SEXP second_counts = PROTECT(Rf_allocVector(INTSXP, count));
        for (R_xlen_t i = 0; i < count; ++i) {
            if (copied[static_cast<std::size_t>(i)] >
                    static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
                retained[static_cast<std::size_t>(i)] >
                    static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
                UNPROTECT(2);
                fail("founder audit count exceeds R integer range");
            }
            INTEGER(first_counts)[i] = static_cast<int>(copied[static_cast<std::size_t>(i)]);
            INTEGER(second_counts)[i] = static_cast<int>(retained[static_cast<std::size_t>(i)]);
        }
        SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(out, 0, first_counts);
        SET_VECTOR_ELT(out, 1, second_counts);
        set_names(out, {"copied_alternative", "retained_alternative"});
        UNPROTECT(3);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("native founder batch materialization: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_make_gamete(
    SEXP destination_pointer, SEXP destination_individual_sexp,
    SEXP parent_h1_pointer, SEXP parent_h2_pointer,
    SEXP parent_individual_sexp, SEXP starting_haplotype_sexp,
    SEXP boundaries) {
    try {
        Packed* destination = require_packed(destination_pointer);
        Packed* h1 = require_packed(parent_h1_pointer);
        Packed* h2 = require_packed(parent_h2_pointer);
        if (TYPEOF(boundaries) != INTSXP) {
            fail("crossover boundaries must be an integer vector");
        }
        std::vector<std::uint64_t> converted(
            static_cast<std::size_t>(XLENGTH(boundaries)));
        for (R_xlen_t i = 0; i < XLENGTH(boundaries); ++i) {
            if (INTEGER(boundaries)[i] == NA_INTEGER ||
                INTEGER(boundaries)[i] < 0) {
                fail("crossover boundaries must be nonnegative");
            }
            converted[static_cast<std::size_t>(i)] =
                static_cast<std::uint64_t>(INTEGER(boundaries)[i]);
        }
        const int starting = scalar_int(starting_haplotype_sexp,
                                        "starting haplotype", 1);
        if (starting > 2) fail("starting haplotype must be 1 or 2");
        destination->value.make_gamete(
            static_cast<std::uint64_t>(scalar_int(
                destination_individual_sexp, "destination individual")),
            h1->value, h2->value,
            static_cast<std::uint64_t>(scalar_int(
                parent_individual_sexp, "parent individual")),
            static_cast<std::uint32_t>(starting - 1),
            converted.empty() ? nullptr : converted.data(),
            static_cast<std::uint64_t>(converted.size()));
        return destination_pointer;
    } catch (const std::exception& ex) {
        Rf_error("native gamete construction: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_decode_genotypes(SEXP h1_pointer,
                                               SEXP h2_pointer) {
    try {
        Packed* h1 = require_packed(h1_pointer);
        Packed* h2 = require_packed(h2_pointer);
        const std::uint64_t individuals = h1->value.individual_count();
        const std::uint64_t markers = h1->value.marker_count();
        if (individuals != h2->value.individual_count() ||
            markers != h2->value.marker_count() ||
            individuals > static_cast<std::uint64_t>(
                                  std::numeric_limits<int>::max()) ||
            markers > static_cast<std::uint64_t>(
                              std::numeric_limits<int>::max())) {
            fail("packed phases have incompatible R dimensions");
        }
        SEXP out = PROTECT(Rf_allocMatrix(
            RAWSXP, static_cast<int>(individuals), static_cast<int>(markers)));
        h1->value.decode_genotypes(h2->value, RAW(out),
                                   static_cast<std::uint64_t>(XLENGTH(out)),
                                   individuals);
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("native bounded genotype decoding: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_bed_sink_create(
    SEXP path_sexp, SEXP individuals_sexp,
    SEXP overwrite_sexp, SEXP buffer_variants_sexp) {
    try {
        const std::string path = scalar_utf8(path_sexp, "BED path");
        const int individuals = scalar_int(individuals_sexp, "individuals", 1);
        const int buffer_variants =
            scalar_int(buffer_variants_sexp, "buffer_variants", 1);
        const bool overwrite = scalar_bool(overwrite_sexp, "overwrite");
        BedSink* sink = new BedSink(
            path, static_cast<std::uint64_t>(individuals), overwrite,
            static_cast<std::uint64_t>(buffer_variants));
        return gsim::native::r::make_owned(sink);
    } catch (const std::exception& ex) {
        Rf_error("native BED sink creation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_bed_sink_append(
    SEXP sink_pointer, SEXP h1_pointer, SEXP h2_pointer) {
    try {
        BedSink* sink = require_bed_sink(sink_pointer);
        Packed* h1 = require_packed(h1_pointer);
        Packed* h2 = require_packed(h2_pointer);
        sink->value.append(h1->value, h2->value);
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("native BED chromosome append: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_bed_sink_finalize(SEXP sink_pointer) {
    try {
        BedSink* sink = require_bed_sink(sink_pointer);
        sink->value.finalize();
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("native BED sink finalization: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_bed_sink_cancel(SEXP sink_pointer) {
    try {
        gsim::native::r::release<BedSink>(sink_pointer, "BED sink");
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("native BED sink cancellation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_bed_sink_info(SEXP sink_pointer) {
    try {
        BedSink* sink = require_bed_sink(sink_pointer);
        SEXP values = PROTECT(Rf_allocVector(REALSXP, 5));
        REAL(values)[0] = static_cast<double>(sink->value.individual_count());
        REAL(values)[1] = static_cast<double>(sink->value.variant_count());
        REAL(values)[2] = static_cast<double>(sink->value.bytes_written());
        REAL(values)[3] = static_cast<double>(sink->value.conversion_buffer_bytes());
        REAL(values)[4] = static_cast<double>(sink->value.lifecycle_object_bytes());
        set_names(values, {"individual_count", "variant_count",
                           "bytes_written", "conversion_buffer_bytes",
                           "lifecycle_object_bytes"});
        const auto state_code = sink->value.state();
        const char* state = state_code == gsim::native::BedSinkState::open
                                ? "open"
                            : state_code == gsim::native::BedSinkState::finalized
                                ? "finalized"
                                : "failed";
        SEXP state_value = PROTECT(Rf_mkString(state));
        Rf_setAttrib(values, Rf_install("state"), state_value);
        UNPROTECT(2);
        return values;
    } catch (const std::exception& ex) {
        Rf_error("native BED sink information: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_sink_create(
    SEXP path_sexp, SEXP individuals_sexp,
    SEXP overwrite_sexp) {
    try {
        const std::string path = scalar_utf8(path_sexp, "HAP path");
        const int individuals = scalar_int(individuals_sexp, "individuals", 1);
        const bool overwrite = scalar_bool(overwrite_sexp, "overwrite");
        HapSink* sink = new HapSink(
            path, static_cast<std::uint64_t>(individuals), overwrite);
        return gsim::native::r::make_owned(sink);
    } catch (const std::exception& ex) {
        Rf_error("native HAP sink creation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_sink_append(
    SEXP sink_pointer, SEXP h1_pointer, SEXP h2_pointer) {
    try {
        HapSink* sink = require_hap_sink(sink_pointer);
        Packed* h1 = require_packed(h1_pointer);
        Packed* h2 = require_packed(h2_pointer);
        sink->value.append(h1->value, h2->value);
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("native HAP append: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_sink_begin(
    SEXP sink_pointer, SEXP marker_count_sexp) {
    try {
        HapSink* sink = require_hap_sink(sink_pointer);
        const int marker_count = scalar_int(marker_count_sexp, "marker_count", 1);
        sink->value.begin_chromosome(static_cast<std::uint64_t>(marker_count));
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("HAP chromosome batch begin: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_sink_write_batch(
    SEXP sink_pointer, SEXP h1_pointer, SEXP h2_pointer,
    SEXP individual_offset_sexp) {
    try {
        HapSink* sink = require_hap_sink(sink_pointer);
        Packed* h1 = require_packed(h1_pointer);
        Packed* h2 = require_packed(h2_pointer);
        const int offset = scalar_int(individual_offset_sexp,
                                      "individual_offset");
        sink->value.write_batch(h1->value, h2->value,
                                static_cast<std::uint64_t>(offset));
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("HAP packed batch write: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_sink_finalize(SEXP sink_pointer) {
    try {
        HapSink* sink = require_hap_sink(sink_pointer);
        sink->value.finalize();
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("native HAP finalization: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_sink_cancel(SEXP sink_pointer) {
    try {
        gsim::native::r::release<HapSink>(sink_pointer, "HAP sink");
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("native HAP cancellation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_sink_info(SEXP sink_pointer) {
    try {
        HapSink* sink = require_hap_sink(sink_pointer);
        SEXP values = PROTECT(Rf_allocVector(REALSXP, 4));
        REAL(values)[0] = exact_r_number(sink->value.individual_count(), "HAP individual count");
        REAL(values)[1] = exact_r_number(sink->value.marker_count(), "HAP marker count");
        REAL(values)[2] = exact_r_number(sink->value.chromosome_count(), "HAP chromosome count");
        REAL(values)[3] = exact_r_number(sink->value.bytes_written(), "HAP byte count");
        set_names(values, {"individual_count", "marker_count",
                           "chromosome_count", "bytes_written"});
        const auto state_code = sink->value.state();
        const char* state = state_code == gsim::native::HapSinkState::open
                                ? "open"
                            : state_code == gsim::native::HapSinkState::finalized
                                ? "finalized"
                                : "failed";
        SEXP state_value = PROTECT(Rf_mkString(state));
        Rf_setAttrib(values, Rf_install("state"), state_value);
        UNPROTECT(2);
        return values;
    } catch (const std::exception& ex) {
        Rf_error("native HAP sink information: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_reader_open(
    SEXP path_sexp) {
    try {
        const std::string path = scalar_utf8(path_sexp, "HAP path");
        HapReader* reader = new HapReader(path);
        return gsim::native::r::make_owned(reader);
    } catch (const std::exception& ex) {
        Rf_error("native HAP reader open: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_reader_close(SEXP reader_pointer) {
    try {
        gsim::native::r::release<HapReader>(reader_pointer, "HAP reader");
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("native HAP reader close: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_reader_info(SEXP reader_pointer) {
    try {
        HapReader* reader = require_hap_reader(reader_pointer);
        const std::uint64_t chromosome_count = reader->value.chromosome_count();
        if (chromosome_count > static_cast<std::uint64_t>(
                std::numeric_limits<int>::max())) fail("too many HAP chromosomes for R");
        const int count = static_cast<int>(chromosome_count);
        SEXP ranges = PROTECT(Rf_allocMatrix(REALSXP, count, 5));
        for (int i = 0; i < count; ++i) {
            const auto& info = reader->value.chromosome(static_cast<std::uint64_t>(i));
            REAL(ranges)[i] = exact_r_number(info.global_start_marker, "HAP marker start");
            REAL(ranges)[i + count] = exact_r_number(info.marker_count, "HAP chromosome markers");
            REAL(ranges)[i + 2 * count] = exact_r_number(info.h1_offset, "HAP H1 offset");
            REAL(ranges)[i + 3 * count] = exact_r_number(info.h2_offset, "HAP H2 offset");
            REAL(ranges)[i + 4 * count] = exact_r_number(info.bytes_per_phase, "HAP phase bytes");
        }
        SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, R_NilValue);
        SEXP columns = PROTECT(Rf_allocVector(STRSXP, 5));
        const char* labels[5] = {"global_start_marker", "marker_count",
                                 "h1_offset", "h2_offset", "bytes_per_phase"};
        for (int i = 0; i < 5; ++i) SET_STRING_ELT(columns, i, Rf_mkChar(labels[i]));
        SET_VECTOR_ELT(dimnames, 1, columns);
        Rf_setAttrib(ranges, R_DimNamesSymbol, dimnames);
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 4));
        SET_VECTOR_ELT(result, 0, Rf_ScalarReal(exact_r_number(reader->value.individual_count(), "HAP individual count")));
        SET_VECTOR_ELT(result, 1, Rf_ScalarReal(exact_r_number(reader->value.marker_count(), "HAP marker count")));
        SET_VECTOR_ELT(result, 2, Rf_ScalarReal(exact_r_number(chromosome_count, "HAP chromosome count")));
        SET_VECTOR_ELT(result, 3, ranges);
        set_names(result, {"individual_count", "marker_count", "chromosome_count", "ranges"});
        UNPROTECT(4);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("native HAP reader information: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_hap_reader_load(
    SEXP reader_pointer, SEXP chromosome_sexp) {
    try {
        HapReader* reader = require_hap_reader(reader_pointer);
        const int chromosome = scalar_int(chromosome_sexp, "chromosome", 1);
        if (static_cast<std::uint64_t>(chromosome) > reader->value.chromosome_count()) {
            fail("HAP chromosome is out of range");
        }
        auto phases = reader->value.load_chromosome(
            static_cast<std::uint64_t>(chromosome - 1));
        SEXP first = PROTECT(make_packed(std::move(phases.first)));
        SEXP second = PROTECT(make_packed(std::move(phases.second)));
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(result, 0, first);
        SET_VECTOR_ELT(result, 1, second);
        set_names(result, {"h1", "h2"});
        UNPROTECT(3);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("native HAP chromosome load: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_packed_bed_read_all(
    SEXP path_sexp, SEXP individuals_sexp,
    SEXP variants_sexp) {
    try {
        const std::string path = scalar_utf8(path_sexp, "BED path");
        const int individuals = scalar_int(individuals_sexp, "individuals", 1);
        const int variants = scalar_int(variants_sexp, "variants", 1);
        BedReader reader(path, static_cast<std::uint64_t>(individuals),
                         static_cast<std::uint64_t>(variants));
        SEXP output = PROTECT(Rf_allocMatrix(INTSXP, individuals, variants));
        std::vector<std::int8_t> record(static_cast<std::size_t>(individuals));
        std::vector<std::uint8_t> bytes(
            static_cast<std::size_t>(reader.bytes_per_variant));
        for (int marker = 0; marker < variants; ++marker) {
            reader.input.read(reinterpret_cast<char*>(bytes.data()),
                              static_cast<std::streamsize>(bytes.size()));
            if (!reader.input) fail("cannot read bounded BED variant");
            for (int individual = 0; individual < individuals; ++individual) {
                const std::uint8_t code =
                    (bytes[static_cast<std::size_t>(individual / 4)] >>
                     (2u * static_cast<unsigned int>(individual % 4))) & 3u;
                static const int dosage[4] = {2, -1, 1, 0};
                INTEGER(output)[individual + individuals * marker] = dosage[code];
            }
        }
        UNPROTECT(1);
        return output;
    } catch (const std::exception& ex) {
        Rf_error("native BED validation decode: %s", ex.what());
    }
    return R_NilValue;
}
