#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace meiosis {

[[noreturn]] void fail(const char* message) {
    throw std::runtime_error(message);
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

private:
    std::uint64_t state_;
};

std::uint64_t fnv1a_utf8(SEXP value, const char* name) {
    if (TYPEOF(value) != STRSXP || XLENGTH(value) != 1 ||
        STRING_ELT(value, 0) == NA_STRING) {
        fail((std::string(name) + " must be one nonmissing string").c_str());
    }
    const unsigned char* text = reinterpret_cast<const unsigned char*>(
        Rf_translateCharUTF8(STRING_ELT(value, 0)));
    std::uint64_t hash = UINT64_C(14695981039346656037);
    while (*text != 0u) {
        hash ^= static_cast<std::uint64_t>(*text++);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

std::uint64_t stream_seed(std::uint64_t seed, SEXP animal, SEXP chromosome,
                          int side) {
    std::uint64_t key = mix64(seed ^ UINT64_C(0x243f6a8885a308d3));
    key ^= mix64(fnv1a_utf8(animal, "animal") ^
                 UINT64_C(0x13198a2e03707344));
    key ^= mix64(fnv1a_utf8(chromosome, "chromosome") ^
                 UINT64_C(0xa4093822299f31d0));
    key ^= side == 1 ? UINT64_C(0x082efa98ec4e6c89)
                     : UINT64_C(0x452821e638d01377);
    return mix64(key);
}

void validate_parent_and_map(SEXP parent_h1, SEXP parent_h2, SEXP positions) {
    if (TYPEOF(parent_h1) != RAWSXP || TYPEOF(parent_h2) != RAWSXP ||
        TYPEOF(positions) != REALSXP || XLENGTH(parent_h1) == 0 ||
        XLENGTH(parent_h1) != XLENGTH(parent_h2) ||
        XLENGTH(parent_h1) != XLENGTH(positions)) {
        fail("parental haplotypes and positions must be nonempty aligned raw/raw/numeric vectors");
    }
    double previous = REAL(positions)[0];
    if (!R_FINITE(previous)) fail("genetic positions must be finite");
    for (R_xlen_t marker = 0; marker < XLENGTH(parent_h1); ++marker) {
        if (RAW(parent_h1)[marker] > 1u || RAW(parent_h2)[marker] > 1u) {
            fail("parental haplotypes must contain only 0/1 alleles");
        }
        const double position = REAL(positions)[marker];
        if (!R_FINITE(position) || (marker > 0 && position < previous)) {
            fail("genetic positions must be finite and nondecreasing");
        }
        previous = position;
    }
}

SEXP materialize(SEXP parent_h1, SEXP parent_h2, SEXP positions,
                 const double* crossovers, R_xlen_t crossover_count,
                 int starting_haplotype) {
    const R_xlen_t marker_count = XLENGTH(parent_h1);
    SEXP gamete = PROTECT(Rf_allocVector(RAWSXP, marker_count));
    int source = starting_haplotype - 1;
    R_xlen_t crossover = 0;
    for (R_xlen_t marker = 0; marker < marker_count; ++marker) {
        while (crossover < crossover_count &&
               crossovers[crossover] <= REAL(positions)[marker]) {
            source ^= 1;
            ++crossover;
        }
        RAW(gamete)[marker] = source == 0 ? RAW(parent_h1)[marker]
                                          : RAW(parent_h2)[marker];
    }
    UNPROTECT(1);
    return gamete;
}

void set_names(SEXP object, const char* first, const char* second,
               const char* third) {
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar(first));
    SET_STRING_ELT(names, 1, Rf_mkChar(second));
    SET_STRING_ELT(names, 2, Rf_mkChar(third));
    Rf_setAttrib(object, R_NamesSymbol, names);
    UNPROTECT(1);
}

}  // namespace meiosis

extern "C" SEXP C_gsim_meiosis_materialize(
    SEXP parent_h1, SEXP parent_h2, SEXP positions, SEXP crossovers,
    SEXP starting_haplotype_sexp) {
    try {
        meiosis::validate_parent_and_map(parent_h1, parent_h2, positions);
        if (TYPEOF(crossovers) != REALSXP) {
            meiosis::fail("crossovers must be a numeric vector");
        }
        if (TYPEOF(starting_haplotype_sexp) != INTSXP ||
            XLENGTH(starting_haplotype_sexp) != 1) {
            meiosis::fail("starting_haplotype must be an integer scalar");
        }
        const int starting_haplotype = INTEGER(starting_haplotype_sexp)[0];
        if (starting_haplotype != 1 && starting_haplotype != 2) {
            meiosis::fail("starting_haplotype must be 1 or 2");
        }
        const double first = REAL(positions)[0];
        const double last = REAL(positions)[XLENGTH(positions) - 1];
        double previous = first;
        for (R_xlen_t i = 0; i < XLENGTH(crossovers); ++i) {
            const double value = REAL(crossovers)[i];
            if (!R_FINITE(value) || value < first || value > last ||
                (i > 0 && value < previous)) {
                meiosis::fail("crossovers must be finite, nondecreasing, and inside the chromosome interval");
            }
            previous = value;
        }
        return meiosis::materialize(parent_h1, parent_h2, positions,
                                    REAL(crossovers), XLENGTH(crossovers),
                                    starting_haplotype);
    } catch (const std::exception& ex) {
        Rf_error("pedigree meiosis materialization: %s", ex.what());
    } catch (...) {
        Rf_error("pedigree meiosis materialization: unknown native error");
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_meiosis_draw(
    SEXP parent_h1, SEXP parent_h2, SEXP positions, SEXP seed_sexp,
    SEXP animal, SEXP chromosome, SEXP side_sexp, SEXP return_audit_sexp) {
    try {
        meiosis::validate_parent_and_map(parent_h1, parent_h2, positions);
        if (TYPEOF(seed_sexp) != REALSXP || XLENGTH(seed_sexp) != 1 ||
            !R_FINITE(REAL(seed_sexp)[0]) || REAL(seed_sexp)[0] < 0.0) {
            meiosis::fail("seed must be one finite nonnegative numeric value");
        }
        if (TYPEOF(side_sexp) != INTSXP || XLENGTH(side_sexp) != 1 ||
            (INTEGER(side_sexp)[0] != 1 && INTEGER(side_sexp)[0] != 2)) {
            meiosis::fail("parental side must be 1 or 2");
        }
        if (TYPEOF(return_audit_sexp) != LGLSXP ||
            XLENGTH(return_audit_sexp) != 1 ||
            LOGICAL(return_audit_sexp)[0] == NA_LOGICAL) {
            meiosis::fail("return_audit must be one logical value");
        }
        const int side = INTEGER(side_sexp)[0];
        meiosis::SplitMix64 rng(meiosis::stream_seed(
            static_cast<std::uint64_t>(REAL(seed_sexp)[0]), animal,
            chromosome, side));
        const double first = REAL(positions)[0];
        const double last = REAL(positions)[XLENGTH(positions) - 1];
        const double length = last - first;

        R_xlen_t crossover_count = 0;
        if (length > 0.0) {
            double arrival = 0.0;
            for (;;) {
                arrival += -std::log(rng.open01());
                if (arrival > length) break;
                if (crossover_count == std::numeric_limits<R_xlen_t>::max()) {
                    meiosis::fail("crossover count exceeds the supported range");
                }
                ++crossover_count;
            }
        }
        SEXP crossovers = PROTECT(Rf_allocVector(REALSXP, crossover_count));
        for (R_xlen_t i = 0; i < crossover_count; ++i) {
            REAL(crossovers)[i] = first + length * rng.open01();
        }
        if (crossover_count > 1) {
            std::sort(REAL(crossovers), REAL(crossovers) + crossover_count);
        }
        const int starting_haplotype = rng.open01() < 0.5 ? 1 : 2;
        SEXP gamete = PROTECT(meiosis::materialize(
            parent_h1, parent_h2, positions, REAL(crossovers),
            crossover_count, starting_haplotype));

        if (LOGICAL(return_audit_sexp)[0] == 0) {
            UNPROTECT(2);
            return gamete;
        }
        SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
        SET_VECTOR_ELT(out, 0, gamete);
        SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(starting_haplotype));
        SET_VECTOR_ELT(out, 2, crossovers);
        meiosis::set_names(out, "gamete", "starting_haplotype", "crossovers");
        UNPROTECT(3);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("pedigree meiosis draw: %s", ex.what());
    } catch (...) {
        Rf_error("pedigree meiosis draw: unknown native error");
    }
    return R_NilValue;
}
