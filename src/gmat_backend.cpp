#include <R.h>
#include <Rinternals.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using status_t = int;
using handle_t = void;

struct WriteInfo {
    std::uint64_t record_count;
    std::uint64_t bytes_written;
    std::uint64_t maximum_record_bytes;
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
    status_t (*variant_write)(const handle_t*, const char*, std::uint32_t,
                              WriteInfo*);
    status_t (*sample_create)(const char* const*, const char* const*,
                              const char* const*, const char* const*,
                              const std::uint32_t*, std::uint64_t, handle_t**);
    status_t (*sample_close)(handle_t*);
    status_t (*sample_write)(const handle_t*, const char*, std::uint32_t,
                             WriteInfo*);
    std::string version;
};

enum class MetadataKind { variant, sample };

struct Metadata {
    Backend* backend;
    handle_t* handle;
    MetadataKind kind;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

template <typename Function>
Function symbol(SEXP symbols, const char* name) {
    SEXP names = Rf_getAttrib(symbols, R_NamesSymbol);
    for (R_xlen_t i = 0; i < XLENGTH(symbols); ++i) {
        if (std::string(CHAR(STRING_ELT(names, i))) == name) {
            SEXP pointer = VECTOR_ELT(symbols, i);
            if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == nullptr) {
                fail(std::string("gmat symbol has no address: ") + name);
            }
            return reinterpret_cast<Function>(R_ExternalPtrAddr(pointer));
        }
    }
    fail(std::string("missing gmat symbol: ") + name);
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

}  // namespace

extern "C" SEXP C_gsim_gmat_backend(SEXP symbols) {
    try {
        if (TYPEOF(symbols) != VECSXP || Rf_getAttrib(symbols, R_NamesSymbol) == R_NilValue) {
            fail("gmat symbols must be a named list");
        }
        Backend* backend = new Backend{};
        try {
            backend->abi_version = symbol<decltype(backend->abi_version)>(symbols, "gmat_abi_version");
            backend->library_version = symbol<decltype(backend->library_version)>(symbols, "gmat_library_version");
            backend->last_error = symbol<decltype(backend->last_error)>(symbols, "gmat_last_error");
            backend->variant_create = symbol<decltype(backend->variant_create)>(symbols, "gmat_variant_metadata_create");
            backend->variant_close = symbol<decltype(backend->variant_close)>(symbols, "gmat_variant_metadata_close");
            backend->variant_write = symbol<decltype(backend->variant_write)>(symbols, "gmat_variant_metadata_write_bim");
            backend->sample_create = symbol<decltype(backend->sample_create)>(symbols, "gmat_sample_metadata_create");
            backend->sample_close = symbol<decltype(backend->sample_close)>(symbols, "gmat_sample_metadata_close");
            backend->sample_write = symbol<decltype(backend->sample_write)>(symbols, "gmat_sample_metadata_write_fam");
            if (backend->abi_version() != 0u) fail("gmat ABI mismatch: gsim requires experimental ABI 0");
            const char* version = backend->library_version();
            if (version == nullptr || version[0] == '\0') fail("gmat returned an empty library version");
            backend->version = version;
        } catch (...) {
            delete backend;
            throw;
        }
        SEXP pointer = PROTECT(R_MakeExternalPtr(backend, R_NilValue, R_NilValue));
        R_RegisterCFinalizerEx(pointer, backend_finalizer, TRUE);
        SEXP version = PROTECT(Rf_mkString(backend->version.c_str()));
        Rf_setAttrib(pointer, Rf_install("gmat_version"), version);
        UNPROTECT(2);
        return pointer;
    } catch (const std::exception& ex) {
        Rf_error("gmat backend: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gmat_variant_create(SEXP backend_pointer, SEXP chromosome,
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

extern "C" SEXP C_gsim_gmat_sample_create(SEXP backend_pointer, SEXP family,
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

extern "C" SEXP C_gsim_gmat_write_bim(SEXP pointer, SEXP path) {
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

extern "C" SEXP C_gsim_gmat_write_fam(SEXP pointer, SEXP path) {
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
