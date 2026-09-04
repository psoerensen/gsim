#include <R.h>
#include <Rinternals.h>

// Native reimplementation of the GPL-3 HAPNEST founder model at commit
// ba52da1a63cf609306ea92540b3d130fa1efd213.  No source was copied; exact
// scientific provenance and intentional interface differences are documented
// in docs/design/hapnest_founder_model.md.

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

[[noreturn]] void fail(const char* message) {
    throw std::runtime_error(message);
}

int scalar_int(SEXP x, const char* name) {
    if (TYPEOF(x) != INTSXP || XLENGTH(x) != 1 || INTEGER(x)[0] == NA_INTEGER) {
        fail((std::string(name) + " must be an integer scalar").c_str());
    }
    return INTEGER(x)[0];
}

double scalar_real(SEXP x, const char* name) {
    if (TYPEOF(x) != REALSXP || XLENGTH(x) != 1 || !R_FINITE(REAL(x)[0])) {
        fail((std::string(name) + " must be a finite numeric scalar").c_str());
    }
    return REAL(x)[0];
}

bool scalar_bool(SEXP x, const char* name) {
    if (TYPEOF(x) != LGLSXP || XLENGTH(x) != 1 || LOGICAL(x)[0] == NA_LOGICAL) {
        fail((std::string(name) + " must be a logical scalar").c_str());
    }
    return LOGICAL(x)[0] != 0;
}

std::uint64_t mix64(std::uint64_t z) noexcept {
    z = (z ^ (z >> 30u)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27u)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31u);
}

class SplitMix64 {
public:
    explicit SplitMix64(std::uint64_t seed) noexcept : state_(seed) {}

    std::uint64_t next() noexcept {
        state_ += UINT64_C(0x9e3779b97f4a7c15);
        return mix64(state_);
    }

    double open01() noexcept {
        const std::uint64_t leading53 = next() >> 11u;
        return (static_cast<double>(leading53) + 0.5) *
               (1.0 / 9007199254740992.0);
    }

    std::uint64_t bounded(std::uint64_t bound) noexcept {
        const std::uint64_t threshold = (UINT64_C(0) - bound) % bound;
        for (;;) {
            const std::uint64_t value = next();
            if (value >= threshold) return value % bound;
        }
    }

private:
    std::uint64_t state_;
};

std::uint64_t stream_seed(std::uint64_t seed, std::uint64_t haplotype,
                          std::uint64_t chromosome_key) noexcept {
    std::uint64_t key = seed ^ UINT64_C(0x6a09e667f3bcc909);
    key ^= (haplotype + 1u) * UINT64_C(0xd2b74407b1ce6e93);
    key ^= chromosome_key * UINT64_C(0xca5a826395121157);
    return mix64(key);
}

std::uint64_t fnv1a_utf8(const std::string& text) noexcept {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned char byte : text) {
        hash ^= static_cast<std::uint64_t>(byte);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

std::vector<std::uint64_t> chromosome_keys(SEXP labels) {
    if (TYPEOF(labels) != STRSXP || XLENGTH(labels) <= 0) {
        fail("chromosome labels must be a nonempty character vector");
    }
    std::vector<std::uint64_t> keys;
    std::vector<std::string> texts;
    keys.reserve(static_cast<std::size_t>(XLENGTH(labels)));
    texts.reserve(static_cast<std::size_t>(XLENGTH(labels)));
    for (R_xlen_t i = 0; i < XLENGTH(labels); ++i) {
        if (STRING_ELT(labels, i) == NA_STRING) {
            fail("chromosome labels must not be missing");
        }
        const std::string text(Rf_translateCharUTF8(STRING_ELT(labels, i)));
        if (text.empty()) fail("chromosome labels must not be empty");
        const std::uint64_t key = fnv1a_utf8(text);
        for (std::size_t previous = 0; previous < keys.size(); ++previous) {
            if (texts[previous] == text) {
                fail("chromosome labels must be unique by exact UTF-8 identity");
            }
            if (keys[previous] == key) {
                fail("distinct chromosome labels collide under FNV-1a");
            }
        }
        texts.push_back(text);
        keys.push_back(key);
    }
    return keys;
}

std::size_t segment_endpoint(const double* position, std::size_t start,
                             std::size_t last, double length) {
    std::size_t endpoint = start;
    double distance = position[endpoint];
    const double objective = distance + length;
    while (distance <= objective && endpoint < last) {
        ++endpoint;
        distance = position[endpoint];
    }
    return endpoint;
}

struct SegmentRecord {
    double individual;
    int phase;
    double haplotype;
    int chromosome_block;
    int start;
    int end;
    int donor_individual;
    int donor_population;
    double coalescent_age;
    double sampled_length;
    double copied_span;
    int copied_alternative;
    int retained_alternative;
};

void set_names(SEXP object, const std::vector<const char*>& names) {
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(names.size())));
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(names.size()); ++i) {
        SET_STRING_ELT(nms, i, Rf_mkChar(names[static_cast<std::size_t>(i)]));
    }
    Rf_setAttrib(object, R_NamesSymbol, nms);
    UNPROTECT(1);
}

SEXP make_segment_columns(const std::vector<SegmentRecord>& records) {
    const R_xlen_t count = static_cast<R_xlen_t>(records.size());
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 13));

    SEXP individual = PROTECT(Rf_allocVector(REALSXP, count));
    SEXP phase = PROTECT(Rf_allocVector(INTSXP, count));
    SEXP haplotype = PROTECT(Rf_allocVector(REALSXP, count));
    SEXP chromosome = PROTECT(Rf_allocVector(INTSXP, count));
    SEXP start = PROTECT(Rf_allocVector(INTSXP, count));
    SEXP end = PROTECT(Rf_allocVector(INTSXP, count));
    SEXP donor = PROTECT(Rf_allocVector(INTSXP, count));
    SEXP population = PROTECT(Rf_allocVector(INTSXP, count));
    SEXP age = PROTECT(Rf_allocVector(REALSXP, count));
    SEXP sampled_length = PROTECT(Rf_allocVector(REALSXP, count));
    SEXP copied_span = PROTECT(Rf_allocVector(REALSXP, count));
    SEXP copied_alternative = PROTECT(Rf_allocVector(INTSXP, count));
    SEXP retained_alternative = PROTECT(Rf_allocVector(INTSXP, count));

    for (R_xlen_t i = 0; i < count; ++i) {
        const SegmentRecord& record = records[static_cast<std::size_t>(i)];
        REAL(individual)[i] = record.individual;
        INTEGER(phase)[i] = record.phase;
        REAL(haplotype)[i] = record.haplotype;
        INTEGER(chromosome)[i] = record.chromosome_block;
        INTEGER(start)[i] = record.start;
        INTEGER(end)[i] = record.end;
        INTEGER(donor)[i] = record.donor_individual;
        INTEGER(population)[i] = record.donor_population;
        REAL(age)[i] = record.coalescent_age;
        REAL(sampled_length)[i] = record.sampled_length;
        REAL(copied_span)[i] = record.copied_span;
        INTEGER(copied_alternative)[i] = record.copied_alternative;
        INTEGER(retained_alternative)[i] = record.retained_alternative;
    }

    SET_VECTOR_ELT(out, 0, individual);
    SET_VECTOR_ELT(out, 1, phase);
    SET_VECTOR_ELT(out, 2, haplotype);
    SET_VECTOR_ELT(out, 3, chromosome);
    SET_VECTOR_ELT(out, 4, start);
    SET_VECTOR_ELT(out, 5, end);
    SET_VECTOR_ELT(out, 6, donor);
    SET_VECTOR_ELT(out, 7, population);
    SET_VECTOR_ELT(out, 8, age);
    SET_VECTOR_ELT(out, 9, sampled_length);
    SET_VECTOR_ELT(out, 10, copied_span);
    SET_VECTOR_ELT(out, 11, copied_alternative);
    SET_VECTOR_ELT(out, 12, retained_alternative);
    set_names(out, {"individual", "phase", "haplotype", "chromosome_block",
                    "start", "end", "donor_individual",
                    "donor_population_code", "coalescent_age",
                    "sampled_length", "copied_genetic_span",
                    "copied_alternative", "retained_alternative"});

    UNPROTECT(14);
    return out;
}

}  // namespace

extern "C" SEXP C_gsim_hapnest_founders(
    SEXP reference_h1, SEXP reference_h2, SEXP donor_codes, SEXP weights,
    SEXP population_n, SEXP population_ne, SEXP population_rho,
    SEXP chromosome_blocks, SEXP chromosome_labels, SEXP positions,
    SEXP mutation_ages,
    SEXP n_individuals_sexp, SEXP seed_sexp, SEXP return_genotypes_sexp,
    SEXP return_segments_sexp, SEXP offset_sexp) {
    try {
        if (TYPEOF(reference_h1) != RAWSXP || !Rf_isMatrix(reference_h1) ||
            TYPEOF(reference_h2) != RAWSXP || !Rf_isMatrix(reference_h2)) {
            fail("H1 and H2 reference haplotypes must be raw matrices");
        }
        SEXP dims = Rf_getAttrib(reference_h1, R_DimSymbol);
        SEXP dims_h2 = Rf_getAttrib(reference_h2, R_DimSymbol);
        if (INTEGER(dims)[0] != INTEGER(dims_h2)[0] ||
            INTEGER(dims)[1] != INTEGER(dims_h2)[1]) {
            fail("H1 and H2 reference haplotypes must have identical dimensions");
        }
        const int donor_count = INTEGER(dims)[0];
        const int marker_count = INTEGER(dims)[1];
        if (donor_count <= 0 || marker_count <= 0) fail("reference matrix is empty");

        if (TYPEOF(donor_codes) != INTSXP || XLENGTH(donor_codes) != donor_count ||
            TYPEOF(weights) != REALSXP || TYPEOF(population_n) != REALSXP ||
            TYPEOF(population_ne) != REALSXP || TYPEOF(population_rho) != REALSXP) {
            fail("invalid population vectors");
        }
        const R_xlen_t population_count = XLENGTH(weights);
        if (population_count <= 0 || XLENGTH(population_n) != population_count ||
            XLENGTH(population_ne) != population_count ||
            XLENGTH(population_rho) != population_count) {
            fail("population vectors have inconsistent lengths");
        }
        if (TYPEOF(chromosome_blocks) != INTSXP ||
            XLENGTH(chromosome_blocks) != marker_count || TYPEOF(positions) != REALSXP ||
            XLENGTH(positions) != marker_count || TYPEOF(mutation_ages) != REALSXP ||
            XLENGTH(mutation_ages) != marker_count) {
            fail("variant metadata vectors are invalid");
        }

        const int n_individuals = scalar_int(n_individuals_sexp, "n");
        const int offset = scalar_int(offset_sexp, "individual_offset");
        const double seed_value = scalar_real(seed_sexp, "seed");
        const bool return_genotypes = scalar_bool(return_genotypes_sexp, "return_genotypes");
        const bool return_segments = scalar_bool(return_segments_sexp, "return_segments");
        if (n_individuals <= 0 || offset < 0 || seed_value < 0.0 ||
            seed_value > 9007199254740991.0) {
            fail("invalid n, offset, or seed");
        }
        const std::uint64_t seed = static_cast<std::uint64_t>(seed_value);

        std::vector<std::vector<int>> donors(static_cast<std::size_t>(population_count));
        for (int donor = 0; donor < donor_count; ++donor) {
            const int code = INTEGER(donor_codes)[donor];
            if (code < 0 || code > population_count) fail("invalid donor population code");
            if (code > 0) donors[static_cast<std::size_t>(code - 1)].push_back(donor);
            for (int marker = 0; marker < marker_count; ++marker) {
                const R_xlen_t index = donor + donor_count * marker;
                if (RAW(reference_h1)[index] > 1u || RAW(reference_h2)[index] > 1u) {
                    fail("reference haplotypes contain a nonbinary allele");
                }
            }
        }

        double weight_sum = 0.0;
        std::vector<double> cumulative(static_cast<std::size_t>(population_count));
        for (R_xlen_t pop = 0; pop < population_count; ++pop) {
            const double w = REAL(weights)[pop];
            const double n_ref = REAL(population_n)[pop];
            const double ne = REAL(population_ne)[pop];
            const double rho = REAL(population_rho)[pop];
            if (!R_FINITE(w) || w <= 0.0 || !R_FINITE(n_ref) || n_ref <= 0.0 ||
                !R_FINITE(ne) || ne <= 0.0 || !R_FINITE(rho) || rho <= 0.0 ||
                donors[static_cast<std::size_t>(pop)].empty()) {
                fail("active population inputs must be finite, positive, and have donors");
            }
            weight_sum += w;
            cumulative[static_cast<std::size_t>(pop)] = weight_sum;
        }
        if (!R_FINITE(weight_sum) || weight_sum <= 0.0) fail("invalid ancestry weights");

        std::vector<std::size_t> block_start;
        std::vector<std::size_t> block_end;
        int previous_block = 0;
        for (int marker = 0; marker < marker_count; ++marker) {
            const int block = INTEGER(chromosome_blocks)[marker];
            const double position = REAL(positions)[marker];
            const double age = REAL(mutation_ages)[marker];
            if (block <= 0 || !R_FINITE(position) || !R_FINITE(age) || age < 0.0) {
                fail("invalid chromosome, position, or mutation age");
            }
            if (marker == 0 || block != previous_block) {
                if (block != previous_block + 1) fail("chromosome blocks must be consecutive");
                if (marker > 0) block_end.push_back(static_cast<std::size_t>(marker - 1));
                block_start.push_back(static_cast<std::size_t>(marker));
                previous_block = block;
            } else if (REAL(positions)[marker] < REAL(positions)[marker - 1]) {
                fail("genetic positions are unsorted within a chromosome");
            }
        }
        block_end.push_back(static_cast<std::size_t>(marker_count - 1));
        const std::vector<std::uint64_t> stable_chromosome_keys =
            chromosome_keys(chromosome_labels);
        if (stable_chromosome_keys.size() != block_start.size()) {
            fail("chromosome label count must equal chromosome block count");
        }

        SEXP h1 = PROTECT(Rf_allocMatrix(RAWSXP, n_individuals, marker_count));
        SEXP h2 = PROTECT(Rf_allocMatrix(RAWSXP, n_individuals, marker_count));
        SEXP genotypes = R_NilValue;
        if (return_genotypes) genotypes = PROTECT(Rf_allocMatrix(RAWSXP, n_individuals, marker_count));
        std::vector<SegmentRecord> records;

        for (int individual = 0; individual < n_individuals; ++individual) {
            const std::uint64_t global_individual =
                static_cast<std::uint64_t>(offset) + static_cast<std::uint64_t>(individual);
            for (int phase = 0; phase < 2; ++phase) {
                Rbyte* output = phase == 0 ? RAW(h1) : RAW(h2);
                const Rbyte* phase_reference =
                    phase == 0 ? RAW(reference_h1) : RAW(reference_h2);
                const std::uint64_t global_haplotype = global_individual * 2u +
                                                       static_cast<std::uint64_t>(phase);
                for (std::size_t block = 0; block < block_start.size(); ++block) {
                    SplitMix64 rng(stream_seed(seed, global_haplotype,
                                              stable_chromosome_keys[block]));
                    std::size_t position = block_start[block];
                    const std::size_t last = block_end[block];
                    while (position <= last) {
                        const double population_draw = rng.open01() * weight_sum;
                        std::size_t pop = 0;
                        while (pop + 1u < cumulative.size() &&
                               population_draw >= cumulative[pop]) {
                            ++pop;
                        }
                        const double gamma_scale = REAL(population_ne)[pop] /
                                                   REAL(population_n)[pop];
                        const double coalescent_age = -gamma_scale *
                            (std::log(rng.open01()) + std::log(rng.open01()));
                        const double exponential_scale =
                            1.0 / (2.0 * coalescent_age * REAL(population_rho)[pop]);
                        const double sampled_length =
                            -exponential_scale * std::log(rng.open01());
                        if (!(coalescent_age > 0.0) || !R_FINITE(coalescent_age) ||
                            !(sampled_length >= 0.0) || ISNAN(sampled_length)) {
                            fail("population parameters produced an invalid random draw");
                        }
                        const std::vector<int>& population_donors = donors[pop];
                        const int donor = population_donors[static_cast<std::size_t>(
                            rng.bounded(static_cast<std::uint64_t>(population_donors.size())))];
                        const std::size_t endpoint = segment_endpoint(
                            REAL(positions), position, last, sampled_length);

                        int copied_alternative = 0;
                        int retained_alternative = 0;
                        for (std::size_t marker = position; marker <= endpoint; ++marker) {
                            const Rbyte allele = phase_reference[
                                donor + donor_count * static_cast<int>(marker)];
                            copied_alternative += allele == 1u ? 1 : 0;
                            const Rbyte retained =
                                allele == 1u && coalescent_age < REAL(mutation_ages)[marker]
                                    ? static_cast<Rbyte>(1u)
                                    : static_cast<Rbyte>(0u);
                            output[individual + n_individuals * static_cast<int>(marker)] = retained;
                            retained_alternative += retained == 1u ? 1 : 0;
                        }

                        if (return_segments) {
                            records.push_back(SegmentRecord{
                                static_cast<double>(global_individual + 1u), phase + 1,
                                static_cast<double>(global_haplotype + 1u),
                                static_cast<int>(block + 1u),
                                static_cast<int>(position + 1u),
                                static_cast<int>(endpoint + 1u), donor + 1,
                                static_cast<int>(pop + 1u), coalescent_age,
                                sampled_length,
                                REAL(positions)[endpoint] - REAL(positions)[position],
                                copied_alternative, retained_alternative});
                        }
                        position = endpoint + 1u;
                    }
                }
            }
        }

        if (return_genotypes) {
            for (R_xlen_t i = 0; i < XLENGTH(h1); ++i) {
                RAW(genotypes)[i] = static_cast<Rbyte>(RAW(h1)[i] + RAW(h2)[i]);
            }
        }

        SEXP segment_columns = R_NilValue;
        if (return_segments) {
            segment_columns = PROTECT(make_segment_columns(records));
        }
        const int output_count = 2 + (return_genotypes ? 1 : 0) +
                                 (return_segments ? 1 : 0);
        SEXP out = PROTECT(Rf_allocVector(VECSXP, output_count));
        SET_VECTOR_ELT(out, 0, h1);
        SET_VECTOR_ELT(out, 1, h2);
        std::vector<const char*> output_names{"h1", "h2"};
        int output_index = 2;
        if (return_genotypes) {
            SET_VECTOR_ELT(out, output_index++, genotypes);
            output_names.push_back("genotypes");
        }
        if (return_segments) {
            SET_VECTOR_ELT(out, output_index, segment_columns);
            output_names.push_back("segments");
        }
        set_names(out, output_names);

        UNPROTECT(3 + (return_genotypes ? 1 : 0) +
                  (return_segments ? 1 : 0));
        return out;
    } catch (const std::exception& ex) {
        Rf_error("HAPNEST founder core: %s", ex.what());
    } catch (...) {
        Rf_error("HAPNEST founder core: unknown native error");
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_hapnest_chromosome_keys(SEXP labels) {
    try {
        const std::vector<std::uint64_t> keys = chromosome_keys(labels);
        if (keys.size() > static_cast<std::size_t>(
                              std::numeric_limits<int>::max())) {
            fail("too many chromosome labels for an R matrix");
        }
        SEXP out = PROTECT(Rf_allocVector(
            RAWSXP, static_cast<R_xlen_t>(keys.size() * 8u)));
        for (std::size_t column = 0; column < keys.size(); ++column) {
            for (unsigned int byte = 0; byte < 8u; ++byte) {
                RAW(out)[static_cast<R_xlen_t>(column * 8u + byte)] =
                    static_cast<Rbyte>((keys[column] >> (byte * 8u)) & 0xffu);
            }
        }
        SEXP dimensions = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(dimensions)[0] = 8;
        INTEGER(dimensions)[1] = static_cast<int>(keys.size());
        Rf_setAttrib(out, R_DimSymbol, dimensions);
        UNPROTECT(2);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("HAPNEST chromosome identity: %s", ex.what());
    } catch (...) {
        Rf_error("HAPNEST chromosome identity: unknown native error");
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_hapnest_segment_endpoint(SEXP positions, SEXP start_sexp,
                                                   SEXP last_sexp, SEXP length_sexp) {
    try {
        if (TYPEOF(positions) != REALSXP || XLENGTH(positions) <= 0) {
            fail("positions must be a numeric vector");
        }
        if (TYPEOF(length_sexp) != REALSXP || XLENGTH(length_sexp) != 1 ||
            ISNAN(REAL(length_sexp)[0])) {
            fail("length must be a nonnegative numeric scalar");
        }
        const int start = scalar_int(start_sexp, "start");
        const int last = scalar_int(last_sexp, "last");
        const double length = REAL(length_sexp)[0];
        if (start < 1 || last < start || last > XLENGTH(positions) ||
            ISNAN(length) || length < 0.0) {
            fail("invalid segment endpoint inputs");
        }
        const std::size_t endpoint = segment_endpoint(
            REAL(positions), static_cast<std::size_t>(start - 1),
            static_cast<std::size_t>(last - 1), length);
        return Rf_ScalarInteger(static_cast<int>(endpoint + 1u));
    } catch (const std::exception& ex) {
        Rf_error("HAPNEST segment endpoint: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_hapnest_copy_segment(SEXP donor, SEXP mutation_ages,
                                               SEXP start_sexp, SEXP end_sexp,
                                               SEXP age_sexp) {
    try {
        if (TYPEOF(donor) != RAWSXP || TYPEOF(mutation_ages) != REALSXP ||
            XLENGTH(donor) != XLENGTH(mutation_ages)) {
            fail("donor and mutation ages must be aligned raw/numeric vectors");
        }
        const int start = scalar_int(start_sexp, "start");
        const int end = scalar_int(end_sexp, "end");
        const double age = scalar_real(age_sexp, "T");
        if (start < 1 || end < start || end > XLENGTH(donor) || age < 0.0) {
            fail("invalid copy-segment inputs");
        }
        SEXP out = PROTECT(Rf_allocVector(RAWSXP, end - start + 1));
        for (int marker = start - 1; marker < end; ++marker) {
            const Rbyte allele = RAW(donor)[marker];
            if (allele > 1u) fail("donor contains a nonbinary allele");
            RAW(out)[marker - start + 1] =
                allele == 1u && age < REAL(mutation_ages)[marker]
                    ? static_cast<Rbyte>(1u)
                    : static_cast<Rbyte>(0u);
        }
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("HAPNEST segment copy: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_hapnest_pair(SEXP h1, SEXP h2) {
    try {
        if (TYPEOF(h1) != RAWSXP || TYPEOF(h2) != RAWSXP ||
            XLENGTH(h1) != XLENGTH(h2)) {
            fail("paired haplotypes must be aligned raw arrays");
        }
        SEXP out = PROTECT(Rf_duplicate(h1));
        for (R_xlen_t i = 0; i < XLENGTH(h1); ++i) {
            if (RAW(h1)[i] > 1u || RAW(h2)[i] > 1u) {
                fail("paired haplotypes contain a nonbinary allele");
            }
            RAW(out)[i] = static_cast<Rbyte>(RAW(h1)[i] + RAW(h2)[i]);
        }
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("HAPNEST haplotype pairing: %s", ex.what());
    }
    return R_NilValue;
}
