#include "bed_reader.h"
#include "native_r.h"
#include <algorithm>
#include <cmath>
#include <climits>
#include <memory>

namespace {
int positive(SEXP x, const char* label) {
    if (TYPEOF(x) != INTSXP || XLENGTH(x) != 1 || INTEGER(x)[0] == NA_INTEGER || INTEGER(x)[0] < 1)
        throw std::runtime_error(std::string(label) + " must be a positive integer");
    return INTEGER(x)[0];
}
std::string path(SEXP x) {
    if (TYPEOF(x) != STRSXP || XLENGTH(x) != 1 || STRING_ELT(x, 0) == NA_STRING)
        throw std::runtime_error("BED path must be a nonmissing string");
    return Rf_translateCharUTF8(STRING_ELT(x, 0));
}
void indices(SEXP x, int upper, const char* label) {
    if (TYPEOF(x) != INTSXP || !XLENGTH(x) || XLENGTH(x) > INT_MAX)
        throw std::runtime_error(std::string(label) + " must be a nonempty integer vector");
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
        if (INTEGER(x)[i] == NA_INTEGER || INTEGER(x)[i] < 1 || INTEGER(x)[i] > upper)
            throw std::runtime_error(std::string(label) + " index out of range");
}
}

extern "C" SEXP C_gsim_bed_read_selected(SEXP file, SEXP samples, SEXP markers,
                                        SEXP rows, SEXP columns) {
    try {
        const int n = positive(samples, "physical samples"), m = positive(markers, "physical markers");
        indices(rows, n, "sample"); indices(columns, m, "marker");
        if (XLENGTH(columns) > 64) throw std::runtime_error("BED decoding is limited to 64 columns per call");
        SEXP out = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(XLENGTH(rows)), static_cast<int>(XLENGTH(columns))));
        gsim::native::BedReader reader(path(file), n, m);
        for (R_xlen_t j = 0; j < XLENGTH(columns); ++j) {
            reader.read_record(INTEGER(columns)[j] - 1);
            for (R_xlen_t i = 0; i < XLENGTH(rows); ++i) {
                const int value = reader.dosage(INTEGER(rows)[i] - 1);
                REAL(out)[j * XLENGTH(rows) + i] = value < 0 ? NA_REAL : value;
            }
        }
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) { Rf_error("native BED selected decode: %s", ex.what()); }
    return R_NilValue;
}

extern "C" SEXP C_gsim_bed_accumulate(SEXP files, SEXP file_index, SEXP columns,
                                     SEXP replacement, SEXP center, SEXP scale, SEXP effects) {
    try {
        if (TYPEOF(files) != VECSXP || !XLENGTH(files) || XLENGTH(files) > INT_MAX)
            throw std::runtime_error("invalid BED file descriptors");
        indices(file_index, static_cast<int>(XLENGTH(files)), "file");
        const R_xlen_t c = XLENGTH(file_index);
        if (TYPEOF(columns) != INTSXP || XLENGTH(columns) != c ||
            TYPEOF(replacement) != REALSXP || XLENGTH(replacement) != c ||
            TYPEOF(center) != REALSXP || XLENGTH(center) != c ||
            TYPEOF(scale) != REALSXP || XLENGTH(scale) != c ||
            TYPEOF(effects) != REALSXP || !Rf_isMatrix(effects))
            throw std::runtime_error("BED effects/transform dimensions do not align");
        SEXP dims = Rf_getAttrib(effects, R_DimSymbol);
        if (INTEGER(dims)[0] != c || INTEGER(dims)[1] < 1)
            throw std::runtime_error("BED effect matrix dimensions do not align");
        const int traits = INTEGER(dims)[1];
        R_xlen_t ns = 0;
        for (R_xlen_t f = 0; f < XLENGTH(files); ++f) {
            SEXP spec = VECTOR_ELT(files, f);
            if (TYPEOF(spec) != VECSXP || XLENGTH(spec) != 4)
                throw std::runtime_error("invalid BED file descriptor");
            (void)path(VECTOR_ELT(spec, 0));
            const int n = positive(VECTOR_ELT(spec, 1), "physical samples");
            (void)positive(VECTOR_ELT(spec, 2), "physical markers");
            indices(VECTOR_ELT(spec, 3), n, "sample");
            if (!f) ns = XLENGTH(VECTOR_ELT(spec, 3));
            if (ns != XLENGTH(VECTOR_ELT(spec, 3))) throw std::runtime_error("BED sample selections differ in length");
        }
        for (R_xlen_t j = 0; j < c; ++j) {
            SEXP spec = VECTOR_ELT(files, INTEGER(file_index)[j] - 1);
            if (INTEGER(columns)[j] == NA_INTEGER || INTEGER(columns)[j] < 1 ||
                INTEGER(columns)[j] > INTEGER(VECTOR_ELT(spec, 2))[0])
                throw std::runtime_error("physical BED marker index out of range");
            if (!R_FINITE(REAL(replacement)[j]) || !R_FINITE(REAL(center)[j]) ||
                !R_FINITE(REAL(scale)[j]) || REAL(scale)[j] <= 0)
                throw std::runtime_error("invalid BED genotype transform");
        }
        for (R_xlen_t i = 0; i < XLENGTH(effects); ++i)
            if (!R_FINITE(REAL(effects)[i])) throw std::runtime_error("nonfinite BED effect");
        SEXP out = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(ns), traits));
        std::fill(REAL(out), REAL(out) + XLENGTH(out), 0.0);
        std::unique_ptr<gsim::native::BedReader> reader;
        int active = -1;
        SEXP rows = R_NilValue;
        for (R_xlen_t j = 0; j < c; ++j) {
            const int f = INTEGER(file_index)[j] - 1;
            if (active != f) {
                reader.reset(); // keep only one physical BED record resident
                SEXP spec = VECTOR_ELT(files, f);
                reader = std::make_unique<gsim::native::BedReader>(path(VECTOR_ELT(spec, 0)),
                    INTEGER(VECTOR_ELT(spec, 1))[0], INTEGER(VECTOR_ELT(spec, 2))[0]);
                rows = VECTOR_ELT(spec, 3); active = f;
            }
            reader->read_record(INTEGER(columns)[j] - 1);
            reader->add_variant(INTEGER(rows), ns, REAL(replacement)[j], REAL(center)[j],
                                REAL(scale)[j], REAL(effects), j, c, traits, REAL(out));
        }
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) { Rf_error("native BED accumulation: %s", ex.what()); }
    return R_NilValue;
}
