#include <R.h>
#include <Rinternals.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "standalone_metadata_api.h"

namespace {

using status_t = int;
using handle_t = void;

struct WriteInfo {
    std::uint64_t record_count;
    std::uint64_t bytes_written;
    std::uint64_t maximum_record_bytes;
};

struct VcfImportInfo {
    std::uint64_t vcf_sample_count;
    std::uint64_t selected_sample_count;
    std::uint64_t total_records_scanned;
    std::uint64_t retained_variants;
    std::uint64_t outside_selected_chromosome;
    std::uint64_t outside_selected_region;
    std::uint64_t indels;
    std::uint64_t multiallelic_records;
    std::uint64_t symbolic_or_breakend_alleles;
    std::uint64_t other_unsupported_alleles;
    std::uint64_t missing_gt;
    std::uint64_t unphased_gt;
    std::uint64_t non_diploid_gt;
    std::uint64_t duplicate_final_ids;
    std::uint64_t maximum_parsing_buffer_bytes;
};

struct Backend {
    std::uint32_t (*abi_version)();
    const char* (*library_version)();
    const char* (*last_error)();
    status_t (*variant_create)(const char* const*, const char* const*,
                               const double*, const std::uint64_t*,
                               const char* const*, const char* const*,
                               std::uint64_t, handle_t**);
    status_t (*variant_close)(handle_t*);
    status_t (*variant_read)(const char*, handle_t**);
    status_t (*variant_count)(const handle_t*, std::uint64_t*);
    status_t (*variant_get)(const handle_t*, std::uint64_t, const char**,
                            const char**, double*, std::uint64_t*,
                            const char**, const char**);
    status_t (*variant_write)(const handle_t*, const char*, std::uint32_t,
                              WriteInfo*);
    status_t (*sample_create)(const char* const*, const char* const*,
                              const char* const*, const char* const*,
                              const std::uint32_t*, std::uint64_t, handle_t**);
    status_t (*sample_close)(handle_t*);
    status_t (*sample_read)(const char*, handle_t**);
    status_t (*sample_count)(const handle_t*, std::uint64_t*);
    status_t (*sample_get)(const handle_t*, std::uint64_t, const char**,
                           const char**, const char**, const char**,
                           std::uint32_t*);
    status_t (*sample_write)(const handle_t*, const char*, std::uint32_t,
                             WriteInfo*);
    status_t (*vcf_open)(const char*, const char* const*, std::uint64_t,
                         const char*, std::uint32_t, std::uint64_t,
                         std::uint64_t, std::uint32_t, handle_t**);
    status_t (*vcf_close)(handle_t*);
    status_t (*vcf_dimensions)(const handle_t*, std::uint64_t*, std::uint64_t*,
                               std::uint64_t*);
    status_t (*vcf_report)(const handle_t*, const char**, VcfImportInfo*);
    status_t (*vcf_sample)(const handle_t*, std::uint64_t, const char**);
    status_t (*vcf_variant)(const handle_t*, std::uint64_t, const char**,
                            const char**, std::uint64_t*, const char**,
                            const char**, std::uint32_t*);
    status_t (*vcf_chromosome)(const handle_t*, std::uint64_t, const char**,
                               std::uint64_t*, std::uint64_t*);
    status_t (*vcf_start)(handle_t*, std::uint64_t);
    status_t (*vcf_next)(handle_t*, std::uint8_t*, std::uint8_t*,
                         std::uint64_t, std::uint64_t*, std::uint32_t*);
    std::string version;
};

enum class MetadataKind { variant, sample };

struct Metadata {
    Backend* backend;
    handle_t* handle;
    MetadataKind kind;
};

struct VcfReader {
    Backend* backend;
    handle_t* handle;
    std::uint64_t samples;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

Backend* require_backend(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == nullptr) {
        fail("invalid gmat backend pointer");
    }
    return static_cast<Backend*>(R_ExternalPtrAddr(pointer));
}

Metadata* require_metadata(SEXP pointer, MetadataKind kind) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == nullptr) {
        fail("invalid gmat metadata pointer");
    }
    Metadata* value = static_cast<Metadata*>(R_ExternalPtrAddr(pointer));
    if (value->handle == nullptr || value->kind != kind) {
        fail("gmat metadata handle is closed or has the wrong kind");
    }
    return value;
}

VcfReader* require_vcf(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == nullptr) {
        fail("invalid gmat VCF reader pointer");
    }
    VcfReader* value = static_cast<VcfReader*>(R_ExternalPtrAddr(pointer));
    if (value->handle == nullptr) fail("gmat VCF reader is closed");
    return value;
}

void check(Backend* backend, status_t status, const char* operation) {
    if (status == 0) return;
    const char* detail = backend->last_error();
    fail(std::string(operation) + ": " +
         (detail != nullptr && detail[0] != '\0' ? detail : "gmat failed"));
}

std::string scalar_utf8(SEXP value, const char* field) {
    if (TYPEOF(value) != STRSXP || XLENGTH(value) != 1 ||
        STRING_ELT(value, 0) == NA_STRING) {
        fail(std::string(field) + " must be one nonmissing string");
    }
    return Rf_translateCharUTF8(STRING_ELT(value, 0));
}

std::vector<std::string> strings(SEXP values, const char* field) {
    if (TYPEOF(values) != STRSXP) {
        fail(std::string(field) + " must be a character vector");
    }
    std::vector<std::string> out(static_cast<std::size_t>(XLENGTH(values)));
    for (R_xlen_t i = 0; i < XLENGTH(values); ++i) {
        if (STRING_ELT(values, i) == NA_STRING) {
            fail(std::string(field) + " must not contain missing values");
        }
        out[static_cast<std::size_t>(i)] =
            Rf_translateCharUTF8(STRING_ELT(values, i));
    }
    return out;
}

std::vector<const char*> pointers(const std::vector<std::string>& values) {
    std::vector<const char*> out(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) out[i] = values[i].c_str();
    return out;
}

void backend_finalizer(SEXP pointer) {
    delete static_cast<Backend*>(R_ExternalPtrAddr(pointer));
    R_ClearExternalPtr(pointer);
}

void metadata_finalizer(SEXP pointer) {
    Metadata* value = static_cast<Metadata*>(R_ExternalPtrAddr(pointer));
    if (value != nullptr) {
        if (value->handle != nullptr) {
            if (value->kind == MetadataKind::variant) {
                (void)value->backend->variant_close(value->handle);
            } else {
                (void)value->backend->sample_close(value->handle);
            }
        }
        delete value;
    }
    R_ClearExternalPtr(pointer);
}

void vcf_finalizer(SEXP pointer) {
    VcfReader* value = static_cast<VcfReader*>(R_ExternalPtrAddr(pointer));
    if (value != nullptr) {
        if (value->handle != nullptr) (void)value->backend->vcf_close(value->handle);
        delete value;
    }
    R_ClearExternalPtr(pointer);
}

SEXP make_info(const WriteInfo& info) {
    SEXP out = PROTECT(Rf_allocVector(REALSXP, 3));
    REAL(out)[0] = static_cast<double>(info.record_count);
    REAL(out)[1] = static_cast<double>(info.bytes_written);
    REAL(out)[2] = static_cast<double>(info.maximum_record_bytes);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("record_count"));
    SET_STRING_ELT(names, 1, Rf_mkChar("bytes_written"));
    SET_STRING_ELT(names, 2, Rf_mkChar("maximum_record_bytes"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

void set_names(SEXP value, const std::vector<const char*>& labels) {
    SEXP names = PROTECT(Rf_allocVector(STRSXP,
        static_cast<R_xlen_t>(labels.size())));
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(labels.size()); ++i) {
        SET_STRING_ELT(names, i, Rf_mkChar(labels[static_cast<std::size_t>(i)]));
    }
    Rf_setAttrib(value, R_NamesSymbol, names);
    UNPROTECT(1);
}

}  // namespace

extern "C" SEXP C_gsim_metadata_backend() {
    try {
        Backend* backend = new Backend{};
        try {
#define GSIM_METADATA_ASSIGN(field, function) \
            backend->field = reinterpret_cast<decltype(backend->field)>(function)
            GSIM_METADATA_ASSIGN(abi_version, native_metadata_abi_version);
            GSIM_METADATA_ASSIGN(library_version, native_metadata_library_version);
            GSIM_METADATA_ASSIGN(last_error, native_metadata_last_error);
            GSIM_METADATA_ASSIGN(variant_create, native_metadata_variant_metadata_create);
            GSIM_METADATA_ASSIGN(variant_close, native_metadata_variant_metadata_close);
            GSIM_METADATA_ASSIGN(variant_read, native_metadata_variant_metadata_read_bim);
            GSIM_METADATA_ASSIGN(variant_count, native_metadata_variant_metadata_count);
            GSIM_METADATA_ASSIGN(variant_get, native_metadata_variant_metadata_get);
            GSIM_METADATA_ASSIGN(variant_write, native_metadata_variant_metadata_write_bim);
            GSIM_METADATA_ASSIGN(sample_create, native_metadata_sample_metadata_create);
            GSIM_METADATA_ASSIGN(sample_close, native_metadata_sample_metadata_close);
            GSIM_METADATA_ASSIGN(sample_read, native_metadata_sample_metadata_read_fam);
            GSIM_METADATA_ASSIGN(sample_count, native_metadata_sample_metadata_count);
            GSIM_METADATA_ASSIGN(sample_get, native_metadata_sample_metadata_get);
            GSIM_METADATA_ASSIGN(sample_write, native_metadata_sample_metadata_write_fam);
            GSIM_METADATA_ASSIGN(vcf_open, native_metadata_phased_vcf_reader_open);
            GSIM_METADATA_ASSIGN(vcf_close, native_metadata_phased_vcf_reader_close);
            GSIM_METADATA_ASSIGN(vcf_dimensions, native_metadata_phased_vcf_reader_dimensions);
            GSIM_METADATA_ASSIGN(vcf_report, native_metadata_phased_vcf_reader_report);
            GSIM_METADATA_ASSIGN(vcf_sample, native_metadata_phased_vcf_reader_sample);
            GSIM_METADATA_ASSIGN(vcf_variant, native_metadata_phased_vcf_reader_variant);
            GSIM_METADATA_ASSIGN(vcf_chromosome, native_metadata_phased_vcf_reader_chromosome);
            GSIM_METADATA_ASSIGN(vcf_start, native_metadata_phased_vcf_reader_start_chromosome);
            GSIM_METADATA_ASSIGN(vcf_next, native_metadata_phased_vcf_reader_next);
#undef GSIM_METADATA_ASSIGN
            backend->version = "gsim-native-1";
        } catch (...) {
            delete backend;
            throw;
        }
        SEXP pointer = PROTECT(R_MakeExternalPtr(backend, R_NilValue, R_NilValue));
        R_RegisterCFinalizerEx(pointer, backend_finalizer, TRUE);
        SEXP version = PROTECT(Rf_mkString(backend->version.c_str()));
        Rf_setAttrib(pointer, Rf_install("native_metadata_version"), version);
        UNPROTECT(2);
        return pointer;
    } catch (const std::exception& ex) {
        Rf_error("gmat backend: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_variant_create(SEXP backend_pointer, SEXP chromosome,
                                             SEXP ids, SEXP cm, SEXP bp,
                                             SEXP alt, SEXP ref) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const auto chromosome_values = strings(chromosome, "chromosome");
        const auto id_values = strings(ids, "variant IDs");
        const auto alt_values = strings(alt, "alternate alleles");
        const auto ref_values = strings(ref, "reference alleles");
        const std::size_t count = chromosome_values.size();
        if (id_values.size() != count || alt_values.size() != count ||
            ref_values.size() != count || TYPEOF(cm) != REALSXP ||
            TYPEOF(bp) != REALSXP || static_cast<std::size_t>(XLENGTH(cm)) != count ||
            static_cast<std::size_t>(XLENGTH(bp)) != count) {
            fail("variant metadata columns must have identical lengths and numeric positions");
        }
        std::vector<std::uint64_t> bp_values(count);
        for (std::size_t i = 0; i < count; ++i) {
            const double value = REAL(bp)[static_cast<R_xlen_t>(i)];
            if (!R_FINITE(value) || value < 1.0 || value > 9007199254740991.0 ||
                std::floor(value) != value) {
                fail("base-pair positions must be exact positive integers not exceeding 2^53-1");
            }
            bp_values[i] = static_cast<std::uint64_t>(value);
        }
        const auto chromosome_ptrs = pointers(chromosome_values);
        const auto id_ptrs = pointers(id_values);
        const auto alt_ptrs = pointers(alt_values);
        const auto ref_ptrs = pointers(ref_values);
        handle_t* handle = nullptr;
        check(backend, backend->variant_create(
              chromosome_ptrs.data(), id_ptrs.data(), REAL(cm), bp_values.data(),
              alt_ptrs.data(), ref_ptrs.data(), static_cast<std::uint64_t>(count),
              &handle), "gmat variant metadata validation");
        Metadata* metadata = new Metadata{backend, handle, MetadataKind::variant};
        SEXP pointer = PROTECT(R_MakeExternalPtr(metadata, R_NilValue, backend_pointer));
        R_RegisterCFinalizerEx(pointer, metadata_finalizer, TRUE);
        UNPROTECT(1);
        return pointer;
    } catch (const std::exception& ex) {
        Rf_error("gmat variant metadata: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_sample_create(SEXP backend_pointer, SEXP family,
                                            SEXP ids, SEXP paternal,
                                            SEXP maternal, SEXP sex) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const auto family_values = strings(family, "family IDs");
        const auto id_values = strings(ids, "individual IDs");
        const auto paternal_values = strings(paternal, "paternal IDs");
        const auto maternal_values = strings(maternal, "maternal IDs");
        const std::size_t count = family_values.size();
        if (id_values.size() != count || paternal_values.size() != count ||
            maternal_values.size() != count || TYPEOF(sex) != INTSXP ||
            static_cast<std::size_t>(XLENGTH(sex)) != count) {
            fail("sample metadata columns must have identical lengths and integer sex");
        }
        std::vector<std::uint32_t> sex_values(count);
        for (std::size_t i = 0; i < count; ++i) {
            const int value = INTEGER(sex)[static_cast<R_xlen_t>(i)];
            if (value == NA_INTEGER || value < 0 || value > 2) {
                fail("sex must be 0 unknown, 1 male, or 2 female");
            }
            sex_values[i] = static_cast<std::uint32_t>(value);
        }
        const auto family_ptrs = pointers(family_values);
        const auto id_ptrs = pointers(id_values);
        const auto paternal_ptrs = pointers(paternal_values);
        const auto maternal_ptrs = pointers(maternal_values);
        handle_t* handle = nullptr;
        check(backend, backend->sample_create(
              family_ptrs.data(), id_ptrs.data(), paternal_ptrs.data(),
              maternal_ptrs.data(), sex_values.data(),
              static_cast<std::uint64_t>(count), &handle),
              "gmat sample metadata validation");
        Metadata* metadata = new Metadata{backend, handle, MetadataKind::sample};
        SEXP pointer = PROTECT(R_MakeExternalPtr(metadata, R_NilValue, backend_pointer));
        R_RegisterCFinalizerEx(pointer, metadata_finalizer, TRUE);
        UNPROTECT(1);
        return pointer;
    } catch (const std::exception& ex) {
        Rf_error("gmat sample metadata: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_write_bim(SEXP pointer, SEXP path) {
    try {
        Metadata* metadata = require_metadata(pointer, MetadataKind::variant);
        const std::string destination = scalar_utf8(path, "BIM path");
        WriteInfo info{};
        check(metadata->backend,
              metadata->backend->variant_write(metadata->handle,
                                                destination.c_str(), 0u, &info),
              "gmat BIM write");
        return make_info(info);
    } catch (const std::exception& ex) {
        Rf_error("gmat BIM write: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_write_fam(SEXP pointer, SEXP path) {
    try {
        Metadata* metadata = require_metadata(pointer, MetadataKind::sample);
        const std::string destination = scalar_utf8(path, "FAM path");
        WriteInfo info{};
        check(metadata->backend,
              metadata->backend->sample_write(metadata->handle,
                                               destination.c_str(), 0u, &info),
              "gmat FAM write");
        return make_info(info);
    } catch (const std::exception& ex) {
        Rf_error("gmat FAM write: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_read_bim(SEXP backend_pointer, SEXP path) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const std::string source = scalar_utf8(path, "BIM path");
        handle_t* handle = nullptr;
        check(backend, backend->variant_read(source.c_str(), &handle),
              "gmat BIM read");
        Metadata* metadata = nullptr;
        try { metadata = new Metadata{backend, handle, MetadataKind::variant}; }
        catch (...) { (void)backend->variant_close(handle); throw; }
        SEXP pointer = PROTECT(R_MakeExternalPtr(metadata, R_NilValue,
                                                 backend_pointer));
        R_RegisterCFinalizerEx(pointer, metadata_finalizer, TRUE);
        std::uint64_t count = 0u;
        check(backend, backend->variant_count(handle, &count),
              "gmat BIM record count");
        if (count > static_cast<std::uint64_t>(R_XLEN_T_MAX)) {
            fail("BIM record count exceeds R vector limits");
        }
        const R_xlen_t length = static_cast<R_xlen_t>(count);
        SEXP chromosome = PROTECT(Rf_allocVector(STRSXP, length));
        SEXP ids = PROTECT(Rf_allocVector(STRSXP, length));
        SEXP cm = PROTECT(Rf_allocVector(REALSXP, length));
        SEXP bp = PROTECT(Rf_allocVector(REALSXP, length));
        SEXP alt = PROTECT(Rf_allocVector(STRSXP, length));
        SEXP ref = PROTECT(Rf_allocVector(STRSXP, length));
        for (std::uint64_t i = 0u; i < count; ++i) {
            const char *chr = nullptr, *id = nullptr, *a1 = nullptr, *a2 = nullptr;
            double genetic = 0.0;
            std::uint64_t physical = 0u;
            check(backend, backend->variant_get(handle, i, &chr, &id, &genetic,
                                                 &physical, &a1, &a2),
                  "gmat BIM record");
            if (physical > 9007199254740991ULL) {
                fail("BIM base-pair position exceeds exact R integer range");
            }
            SET_STRING_ELT(chromosome, static_cast<R_xlen_t>(i), Rf_mkCharCE(chr, CE_UTF8));
            SET_STRING_ELT(ids, static_cast<R_xlen_t>(i), Rf_mkCharCE(id, CE_UTF8));
            REAL(cm)[static_cast<R_xlen_t>(i)] = genetic;
            REAL(bp)[static_cast<R_xlen_t>(i)] = static_cast<double>(physical);
            SET_STRING_ELT(alt, static_cast<R_xlen_t>(i), Rf_mkCharCE(a1, CE_UTF8));
            SET_STRING_ELT(ref, static_cast<R_xlen_t>(i), Rf_mkCharCE(a2, CE_UTF8));
        }
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 7));
        SET_VECTOR_ELT(result, 0, pointer); SET_VECTOR_ELT(result, 1, chromosome);
        SET_VECTOR_ELT(result, 2, ids); SET_VECTOR_ELT(result, 3, cm);
        SET_VECTOR_ELT(result, 4, bp); SET_VECTOR_ELT(result, 5, alt);
        SET_VECTOR_ELT(result, 6, ref);
        set_names(result, {"pointer", "chromosome", "variant_id",
                           "genetic_position_cm", "base_pair_position",
                           "alt", "ref"});
        UNPROTECT(8);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("gmat BIM read: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_read_fam(SEXP backend_pointer, SEXP path) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const std::string source = scalar_utf8(path, "FAM path");
        handle_t* handle = nullptr;
        check(backend, backend->sample_read(source.c_str(), &handle),
              "gmat FAM read");
        Metadata* metadata = nullptr;
        try { metadata = new Metadata{backend, handle, MetadataKind::sample}; }
        catch (...) { (void)backend->sample_close(handle); throw; }
        SEXP pointer = PROTECT(R_MakeExternalPtr(metadata, R_NilValue,
                                                 backend_pointer));
        R_RegisterCFinalizerEx(pointer, metadata_finalizer, TRUE);
        std::uint64_t count = 0u;
        check(backend, backend->sample_count(handle, &count),
              "gmat FAM record count");
        if (count > static_cast<std::uint64_t>(R_XLEN_T_MAX)) {
            fail("FAM record count exceeds R vector limits");
        }
        const R_xlen_t length = static_cast<R_xlen_t>(count);
        SEXP family = PROTECT(Rf_allocVector(STRSXP, length));
        SEXP ids = PROTECT(Rf_allocVector(STRSXP, length));
        SEXP paternal = PROTECT(Rf_allocVector(STRSXP, length));
        SEXP maternal = PROTECT(Rf_allocVector(STRSXP, length));
        SEXP sex = PROTECT(Rf_allocVector(INTSXP, length));
        for (std::uint64_t i = 0u; i < count; ++i) {
            const char *fid = nullptr, *id = nullptr, *sire = nullptr, *dam = nullptr;
            std::uint32_t sex_value = 0u;
            check(backend, backend->sample_get(handle, i, &fid, &id, &sire,
                                                &dam, &sex_value),
                  "gmat FAM record");
            SET_STRING_ELT(family, static_cast<R_xlen_t>(i), Rf_mkCharCE(fid, CE_UTF8));
            SET_STRING_ELT(ids, static_cast<R_xlen_t>(i), Rf_mkCharCE(id, CE_UTF8));
            SET_STRING_ELT(paternal, static_cast<R_xlen_t>(i), Rf_mkCharCE(sire, CE_UTF8));
            SET_STRING_ELT(maternal, static_cast<R_xlen_t>(i), Rf_mkCharCE(dam, CE_UTF8));
            INTEGER(sex)[static_cast<R_xlen_t>(i)] = static_cast<int>(sex_value);
        }
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 6));
        SET_VECTOR_ELT(result, 0, pointer); SET_VECTOR_ELT(result, 1, family);
        SET_VECTOR_ELT(result, 2, ids); SET_VECTOR_ELT(result, 3, paternal);
        SET_VECTOR_ELT(result, 4, maternal); SET_VECTOR_ELT(result, 5, sex);
        set_names(result, {"pointer", "family_id", "individual_id",
                           "paternal_id", "maternal_id", "sex"});
        UNPROTECT(7);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("gmat FAM read: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_vcf_open(
    SEXP backend_pointer, SEXP path, SEXP selected_samples,
    SEXP selected_chromosome, SEXP region, SEXP unsupported) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const std::string source = scalar_utf8(path, "VCF path");
        std::vector<std::string> requested_samples;
        if (selected_samples != R_NilValue) {
            requested_samples = strings(selected_samples, "selected samples");
            if (requested_samples.empty()) fail("selected samples must not be empty");
        }
        const auto requested_pointers = pointers(requested_samples);
        std::string chromosome_value;
        if (selected_chromosome != R_NilValue) {
            chromosome_value = scalar_utf8(selected_chromosome, "selected chromosome");
            if (chromosome_value.empty()) fail("selected chromosome must not be empty");
        }
        bool has_region = region != R_NilValue;
        std::uint64_t region_values[2] = {0u, 0u};
        if (has_region) {
            if ((TYPEOF(region) != REALSXP && TYPEOF(region) != INTSXP) ||
                XLENGTH(region) != 2) {
                fail("region must be a numeric vector of length two");
            }
            for (R_xlen_t i = 0; i < 2; ++i) {
                const double value = TYPEOF(region) == REALSXP ? REAL(region)[i] :
                    (INTEGER(region)[i] == NA_INTEGER ? NA_REAL :
                                                       static_cast<double>(INTEGER(region)[i]));
                if (!std::isfinite(value) || value <= 0.0 || value != std::floor(value) ||
                    value > 9007199254740991.0) {
                    fail("region bounds must be exact positive integers");
                }
                region_values[i] = static_cast<std::uint64_t>(value);
            }
            if (chromosome_value.empty() || region_values[0] > region_values[1]) {
                fail("region requires a chromosome and start <= end");
            }
        }
        const std::string unsupported_value = scalar_utf8(unsupported, "unsupported");
        if (unsupported_value != "skip" && unsupported_value != "error") {
            fail("unsupported must be 'skip' or 'error'");
        }
        handle_t* handle = nullptr;
        check(backend, backend->vcf_open(
                  source.c_str(),
                  requested_pointers.empty() ? nullptr : requested_pointers.data(),
                  static_cast<std::uint64_t>(requested_pointers.size()),
                  chromosome_value.empty() ? nullptr : chromosome_value.c_str(),
                  has_region ? 1u : 0u, region_values[0], region_values[1],
                  unsupported_value == "skip" ? 1u : 0u, &handle),
              "native VCF open");
        std::uint64_t sample_count = 0, variant_count = 0, chromosome_count = 0;
        check(backend, backend->vcf_dimensions(handle, &sample_count, &variant_count,
                                               &chromosome_count),
              "gmat VCF dimensions");
        const char* input_type = nullptr;
        VcfImportInfo import_info{};
        check(backend, backend->vcf_report(handle, &input_type, &import_info),
              "native VCF report");
        if (sample_count > static_cast<std::uint64_t>(R_XLEN_T_MAX) ||
            variant_count > static_cast<std::uint64_t>(R_XLEN_T_MAX) ||
            chromosome_count > static_cast<std::uint64_t>(R_XLEN_T_MAX)) {
            (void)backend->vcf_close(handle);
            fail("VCF dimensions exceed R vector limits");
        }
        VcfReader* reader = nullptr;
        try { reader = new VcfReader{backend, handle, sample_count}; }
        catch (...) { (void)backend->vcf_close(handle); throw; }
        SEXP pointer = PROTECT(R_MakeExternalPtr(reader, R_NilValue, backend_pointer));
        R_RegisterCFinalizerEx(pointer, vcf_finalizer, TRUE);
        SEXP samples = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(sample_count)));
        for (std::uint64_t i = 0; i < sample_count; ++i) {
            const char* value = nullptr;
            check(backend, backend->vcf_sample(handle, i, &value), "gmat VCF sample");
            SET_STRING_ELT(samples, static_cast<R_xlen_t>(i), Rf_mkCharCE(value, CE_UTF8));
        }
        SEXP chromosome = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(variant_count)));
        SEXP ids = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(variant_count)));
        SEXP bp = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(variant_count)));
        SEXP ref = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(variant_count)));
        SEXP alt = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(variant_count)));
        SEXP generated = PROTECT(Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(variant_count)));
        for (std::uint64_t i = 0; i < variant_count; ++i) {
            const char *chr = nullptr, *id = nullptr, *reference = nullptr, *alternate = nullptr;
            std::uint64_t physical = 0;
            std::uint32_t made = 0;
            check(backend, backend->vcf_variant(handle, i, &chr, &id, &physical,
                                                 &reference, &alternate, &made),
                  "gmat VCF variant");
            if (physical > 9007199254740991ULL) fail("VCF POS exceeds exact R numeric range");
            SET_STRING_ELT(chromosome, static_cast<R_xlen_t>(i), Rf_mkCharCE(chr, CE_UTF8));
            SET_STRING_ELT(ids, static_cast<R_xlen_t>(i), Rf_mkCharCE(id, CE_UTF8));
            REAL(bp)[static_cast<R_xlen_t>(i)] = static_cast<double>(physical);
            SET_STRING_ELT(ref, static_cast<R_xlen_t>(i), Rf_mkCharCE(reference, CE_UTF8));
            SET_STRING_ELT(alt, static_cast<R_xlen_t>(i), Rf_mkCharCE(alternate, CE_UTF8));
            LOGICAL(generated)[static_cast<R_xlen_t>(i)] = made ? TRUE : FALSE;
        }
        SEXP block_label = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(chromosome_count)));
        SEXP block_first = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(chromosome_count)));
        SEXP block_count = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(chromosome_count)));
        for (std::uint64_t i = 0; i < chromosome_count; ++i) {
            const char* label = nullptr;
            std::uint64_t first = 0, count = 0;
            check(backend, backend->vcf_chromosome(handle, i, &label, &first, &count),
                  "gmat VCF chromosome");
            SET_STRING_ELT(block_label, static_cast<R_xlen_t>(i), Rf_mkCharCE(label, CE_UTF8));
            REAL(block_first)[static_cast<R_xlen_t>(i)] = static_cast<double>(first + 1u);
            REAL(block_count)[static_cast<R_xlen_t>(i)] = static_cast<double>(count);
        }
        SEXP variants = PROTECT(Rf_allocVector(VECSXP, 6));
        SET_VECTOR_ELT(variants, 0, chromosome); SET_VECTOR_ELT(variants, 1, ids);
        SET_VECTOR_ELT(variants, 2, bp); SET_VECTOR_ELT(variants, 3, ref);
        SET_VECTOR_ELT(variants, 4, alt); SET_VECTOR_ELT(variants, 5, generated);
        set_names(variants, {"chromosome", "variant_id", "base_pair_position",
                             "ref", "alt", "generated_id"});
        SEXP blocks = PROTECT(Rf_allocVector(VECSXP, 3));
        SET_VECTOR_ELT(blocks, 0, block_label); SET_VECTOR_ELT(blocks, 1, block_first);
        SET_VECTOR_ELT(blocks, 2, block_count);
        set_names(blocks, {"chromosome", "first_variant", "variant_count"});
        SEXP report = PROTECT(Rf_allocVector(VECSXP, 16));
        SET_VECTOR_ELT(report, 0, Rf_mkString(input_type));
        const std::uint64_t report_values[] = {
            import_info.vcf_sample_count, import_info.selected_sample_count,
            import_info.total_records_scanned, import_info.retained_variants,
            import_info.outside_selected_chromosome, import_info.outside_selected_region,
            import_info.indels, import_info.multiallelic_records,
            import_info.symbolic_or_breakend_alleles,
            import_info.other_unsupported_alleles, import_info.missing_gt,
            import_info.unphased_gt, import_info.non_diploid_gt,
            import_info.duplicate_final_ids, import_info.maximum_parsing_buffer_bytes};
        for (R_xlen_t i = 0; i < 15; ++i) {
            SET_VECTOR_ELT(report, i + 1, Rf_ScalarReal(
                static_cast<double>(report_values[static_cast<std::size_t>(i)])));
        }
        set_names(report, {"input_type", "vcf_sample_count", "selected_sample_count",
                           "total_records_scanned", "retained_variants",
                           "outside_selected_chromosome", "outside_selected_region",
                           "indels", "multiallelic_records",
                           "symbolic_or_breakend_alleles", "other_unsupported_alleles",
                           "missing_gt", "unphased_gt", "non_diploid_gt",
                           "duplicate_final_ids", "maximum_parsing_buffer_bytes"});
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 5));
        SET_VECTOR_ELT(result, 0, pointer); SET_VECTOR_ELT(result, 1, samples);
        SET_VECTOR_ELT(result, 2, variants); SET_VECTOR_ELT(result, 3, blocks);
        SET_VECTOR_ELT(result, 4, report);
        set_names(result, {"pointer", "samples", "variants", "chromosomes", "report"});
        UNPROTECT(15);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("gmat VCF open: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_vcf_start(SEXP pointer, SEXP chromosome_index) {
    try {
        VcfReader* reader = require_vcf(pointer);
        if (TYPEOF(chromosome_index) != INTSXP || XLENGTH(chromosome_index) != 1 ||
            INTEGER(chromosome_index)[0] == NA_INTEGER || INTEGER(chromosome_index)[0] < 1) {
            fail("chromosome index must be a positive integer");
        }
        check(reader->backend,
              reader->backend->vcf_start(reader->handle,
                  static_cast<std::uint64_t>(INTEGER(chromosome_index)[0] - 1)),
              "gmat VCF chromosome start");
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("gmat VCF start: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_vcf_next(SEXP pointer) {
    try {
        VcfReader* reader = require_vcf(pointer);
        if (reader->samples > static_cast<std::uint64_t>(R_XLEN_T_MAX)) {
            fail("VCF sample count exceeds R vector limits");
        }
        const R_xlen_t count = static_cast<R_xlen_t>(reader->samples);
        SEXP h1 = PROTECT(Rf_allocVector(RAWSXP, count));
        SEXP h2 = PROTECT(Rf_allocVector(RAWSXP, count));
        std::uint64_t variant = 0;
        std::uint32_t has_record = 0;
        check(reader->backend,
              reader->backend->vcf_next(reader->handle, RAW(h1), RAW(h2),
                                         reader->samples, &variant, &has_record),
              "gmat VCF record");
        if (!has_record) {
            UNPROTECT(2);
            return R_NilValue;
        }
        if (variant >= 9007199254740991ULL) fail("VCF variant index exceeds exact R numeric range");
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
        SET_VECTOR_ELT(result, 0, Rf_ScalarReal(static_cast<double>(variant + 1u)));
        SET_VECTOR_ELT(result, 1, h1); SET_VECTOR_ELT(result, 2, h2);
        set_names(result, {"variant_index", "h1", "h2"});
        UNPROTECT(3);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("gmat VCF next: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_vcf_close(SEXP pointer) {
    try {
        VcfReader* reader = require_vcf(pointer);
        check(reader->backend, reader->backend->vcf_close(reader->handle),
              "gmat VCF close");
        reader->handle = nullptr;
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("gmat VCF close: %s", ex.what());
    }
    return R_NilValue;
}
