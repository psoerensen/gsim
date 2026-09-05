#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using status_t = int;
using handle_t = void;

struct Backend {
    std::uint32_t (*abi_version)();
    const char* (*library_version)();
    const char* (*last_error)();
    status_t (*create_zero)(std::uint64_t, std::uint64_t, handle_t**);
    status_t (*create_values)(std::uint64_t, std::uint64_t,
                              const std::uint8_t*, std::uint64_t,
                              std::uint64_t, handle_t**);
    status_t (*close)(handle_t*);
    status_t (*individual_count)(const handle_t*, std::uint64_t*);
    status_t (*marker_count)(const handle_t*, std::uint64_t*);
    status_t (*words_per_marker)(const handle_t*, std::uint64_t*);
    status_t (*storage_bytes)(const handle_t*, std::uint64_t*);
    status_t (*word)(const handle_t*, std::uint64_t, std::uint64_t,
                     std::uint64_t*);
    status_t (*allele)(const handle_t*, std::uint64_t, std::uint64_t,
                       std::uint8_t*);
    status_t (*set_allele)(handle_t*, std::uint64_t, std::uint64_t,
                           std::uint8_t);
    status_t (*unpack)(const handle_t*, std::uint8_t*, std::uint64_t,
                       std::uint64_t);
    status_t (*copy_interval)(handle_t*, std::uint64_t, const handle_t*,
                              std::uint64_t, std::uint64_t, std::uint64_t);
    status_t (*copy_filtered)(handle_t*, std::uint64_t, const handle_t*,
                              std::uint64_t, std::uint64_t, std::uint64_t,
                              double, const double*, std::uint64_t);
    status_t (*make_gamete)(handle_t*, std::uint64_t, const handle_t*,
                            const handle_t*, std::uint64_t, std::uint32_t,
                            const std::uint64_t*, std::uint64_t);
    status_t (*decode_genotypes)(const handle_t*, const handle_t*,
                                 std::uint8_t*, std::uint64_t,
                                 std::uint64_t);
    std::string version;
};

struct Packed {
    Backend* backend;
    handle_t* handle;
    std::uint64_t individuals;
    std::uint64_t markers;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

SEXP named_element(SEXP values, const char* name) {
    if (TYPEOF(values) != VECSXP) fail("gbits symbol table must be a list");
    SEXP names = Rf_getAttrib(values, R_NamesSymbol);
    if (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(values)) {
        fail("gbits symbol table must be named");
    }
    for (R_xlen_t i = 0; i < XLENGTH(values); ++i) {
        if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            return VECTOR_ELT(values, i);
        }
    }
    fail(std::string("gbits library is missing required symbol '") + name + "'");
}

template <typename Function>
Function symbol(SEXP values, const char* name) {
    SEXP address = named_element(values, name);
    if (TYPEOF(address) != EXTPTRSXP || R_ExternalPtrAddr(address) == nullptr) {
        fail(std::string("gbits symbol '") + name + "' has no callable address");
    }
    return reinterpret_cast<Function>(R_ExternalPtrAddr(address));
}

void backend_finalizer(SEXP pointer) {
    Backend* backend = static_cast<Backend*>(R_ExternalPtrAddr(pointer));
    delete backend;
    R_ClearExternalPtr(pointer);
}

void packed_finalizer(SEXP pointer) {
    Packed* packed = static_cast<Packed*>(R_ExternalPtrAddr(pointer));
    if (packed != nullptr) {
        if (packed->handle != nullptr && packed->backend != nullptr) {
            (void)packed->backend->close(packed->handle);
        }
        delete packed;
    }
    R_ClearExternalPtr(pointer);
}

Backend* require_backend(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP) fail("gbits backend is invalid");
    Backend* backend = static_cast<Backend*>(R_ExternalPtrAddr(pointer));
    if (backend == nullptr) fail("gbits backend has been released");
    return backend;
}

Packed* require_packed(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP) fail("packed haplotypes are invalid");
    Packed* packed = static_cast<Packed*>(R_ExternalPtrAddr(pointer));
    if (packed == nullptr || packed->handle == nullptr) {
        fail("packed haplotypes have been released");
    }
    return packed;
}

void require_same_backend(const Packed* first, const Packed* second) {
    if (first->backend != second->backend) {
        fail("packed haplotypes originate from different gbits backends");
    }
}

void check(Backend* backend, status_t status, const char* operation) {
    if (status == 0) return;
    const char* detail = backend->last_error();
    fail(std::string(operation) + " failed" +
         (detail != nullptr && detail[0] != '\0'
              ? std::string(": ") + detail
              : std::string()));
}

int scalar_int(SEXP value, const char* name, int lower = 0) {
    if (TYPEOF(value) != INTSXP || XLENGTH(value) != 1 ||
        INTEGER(value)[0] == NA_INTEGER || INTEGER(value)[0] < lower) {
        fail(std::string(name) + " must be an integer scalar");
    }
    return INTEGER(value)[0];
}

SEXP make_packed(SEXP backend_pointer, handle_t* handle,
                 std::uint64_t individuals, std::uint64_t markers) {
    Backend* backend = require_backend(backend_pointer);
    Packed* packed = nullptr;
    try {
        packed = new Packed{backend, handle, individuals, markers};
    } catch (...) {
        (void)backend->close(handle);
        throw;
    }
    SEXP pointer = PROTECT(R_MakeExternalPtr(packed, R_NilValue,
                                             backend_pointer));
    R_RegisterCFinalizerEx(pointer, packed_finalizer, TRUE);
    UNPROTECT(1);
    return pointer;
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

} // namespace

extern "C" SEXP C_gsim_gbits_backend(SEXP symbols) {
    try {
        Backend* backend = new Backend{};
        try {
            backend->abi_version = symbol<decltype(backend->abi_version)>(
                symbols, "gbits_abi_version");
            backend->library_version = symbol<decltype(backend->library_version)>(
                symbols, "gbits_library_version");
            backend->last_error = symbol<decltype(backend->last_error)>(
                symbols, "gbits_last_error");
            backend->create_zero = symbol<decltype(backend->create_zero)>(
                symbols, "gbits_phased_haplotype_create_zero");
            backend->create_values = symbol<decltype(backend->create_values)>(
                symbols, "gbits_phased_haplotype_create_from_values");
            backend->close = symbol<decltype(backend->close)>(
                symbols, "gbits_phased_haplotype_close");
            backend->individual_count = symbol<decltype(backend->individual_count)>(
                symbols, "gbits_phased_haplotype_individual_count");
            backend->marker_count = symbol<decltype(backend->marker_count)>(
                symbols, "gbits_phased_haplotype_marker_count");
            backend->words_per_marker = symbol<decltype(backend->words_per_marker)>(
                symbols, "gbits_phased_haplotype_words_per_marker");
            backend->storage_bytes = symbol<decltype(backend->storage_bytes)>(
                symbols, "gbits_phased_haplotype_storage_bytes");
            backend->word = symbol<decltype(backend->word)>(
                symbols, "gbits_phased_haplotype_word");
            backend->allele = symbol<decltype(backend->allele)>(
                symbols, "gbits_phased_haplotype_allele");
            backend->set_allele = symbol<decltype(backend->set_allele)>(
                symbols, "gbits_phased_haplotype_set_allele");
            backend->unpack = symbol<decltype(backend->unpack)>(
                symbols, "gbits_phased_haplotype_unpack");
            backend->copy_interval = symbol<decltype(backend->copy_interval)>(
                symbols, "gbits_phased_haplotype_copy_interval");
            backend->copy_filtered = symbol<decltype(backend->copy_filtered)>(
                symbols, "gbits_phased_haplotype_copy_filtered_segment");
            backend->make_gamete = symbol<decltype(backend->make_gamete)>(
                symbols, "gbits_phased_haplotype_make_gamete");
            backend->decode_genotypes = symbol<decltype(backend->decode_genotypes)>(
                symbols, "gbits_phased_haplotype_decode_genotypes");
            const std::uint32_t abi = backend->abi_version();
            if (abi != 4u) {
                fail("gbits ABI mismatch: gsim requires ABI 4");
            }
            const char* version = backend->library_version();
            if (version == nullptr || version[0] == '\0') {
                fail("gbits returned an empty library version");
            }
            backend->version = version;
        } catch (...) {
            delete backend;
            throw;
        }
        SEXP pointer = PROTECT(R_MakeExternalPtr(backend, R_NilValue,
                                                 R_NilValue));
        R_RegisterCFinalizerEx(pointer, backend_finalizer, TRUE);
        SEXP version = PROTECT(Rf_mkString(backend->version.c_str()));
        Rf_setAttrib(pointer, Rf_install("gbits_version"), version);
        UNPROTECT(2);
        return pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits backend: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_pack(SEXP backend_pointer, SEXP values) {
    try {
        Backend* backend = require_backend(backend_pointer);
        if (TYPEOF(values) != RAWSXP || !Rf_isMatrix(values)) {
            fail("packed input must be a raw matrix");
        }
        SEXP dimensions = Rf_getAttrib(values, R_DimSymbol);
        const int individuals = INTEGER(dimensions)[0];
        const int markers = INTEGER(dimensions)[1];
        if (individuals <= 0 || markers <= 0) fail("packed input is empty");
        handle_t* handle = nullptr;
        check(backend, backend->create_values(
                           static_cast<std::uint64_t>(individuals),
                           static_cast<std::uint64_t>(markers), RAW(values),
                           static_cast<std::uint64_t>(XLENGTH(values)),
                           static_cast<std::uint64_t>(individuals), &handle),
              "gbits pack");
        return make_packed(backend_pointer, handle,
                           static_cast<std::uint64_t>(individuals),
                           static_cast<std::uint64_t>(markers));
    } catch (const std::exception& ex) {
        Rf_error("gbits pack: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_zero(SEXP backend_pointer, SEXP individuals_sexp,
                                   SEXP markers_sexp) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const int individuals = scalar_int(individuals_sexp, "individuals", 1);
        const int markers = scalar_int(markers_sexp, "markers", 1);
        handle_t* handle = nullptr;
        check(backend, backend->create_zero(
                           static_cast<std::uint64_t>(individuals),
                           static_cast<std::uint64_t>(markers), &handle),
              "gbits zero allocation");
        return make_packed(backend_pointer, handle,
                           static_cast<std::uint64_t>(individuals),
                           static_cast<std::uint64_t>(markers));
    } catch (const std::exception& ex) {
        Rf_error("gbits allocation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_unpack(SEXP pointer) {
    try {
        Packed* packed = require_packed(pointer);
        if (packed->individuals > static_cast<std::uint64_t>(
                                      std::numeric_limits<int>::max()) ||
            packed->markers > static_cast<std::uint64_t>(
                                  std::numeric_limits<int>::max())) {
            fail("packed dimensions exceed R matrix limits");
        }
        SEXP out = PROTECT(Rf_allocMatrix(
            RAWSXP, static_cast<int>(packed->individuals),
            static_cast<int>(packed->markers)));
        check(packed->backend,
              packed->backend->unpack(
                  packed->handle, RAW(out),
                  static_cast<std::uint64_t>(XLENGTH(out)),
                  packed->individuals),
              "gbits unpack");
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("gbits unpack: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_info(SEXP pointer) {
    try {
        Packed* packed = require_packed(pointer);
        std::uint64_t words = 0u;
        std::uint64_t bytes = 0u;
        check(packed->backend,
              packed->backend->words_per_marker(packed->handle, &words),
              "gbits word-count query");
        check(packed->backend,
              packed->backend->storage_bytes(packed->handle, &bytes),
              "gbits storage-byte query");
        SEXP out = PROTECT(Rf_allocVector(REALSXP, 4));
        REAL(out)[0] = static_cast<double>(packed->individuals);
        REAL(out)[1] = static_cast<double>(packed->markers);
        REAL(out)[2] = static_cast<double>(words);
        REAL(out)[3] = static_cast<double>(bytes);
        set_names(out, {"individuals", "markers", "words_per_marker",
                        "storage_bytes"});
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("gbits info: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_word(SEXP pointer, SEXP marker_sexp,
                                   SEXP word_sexp) {
    try {
        Packed* packed = require_packed(pointer);
        const int marker = scalar_int(marker_sexp, "marker");
        const int word_index = scalar_int(word_sexp, "word_index");
        std::uint64_t value = 0u;
        check(packed->backend,
              packed->backend->word(packed->handle,
                                    static_cast<std::uint64_t>(marker),
                                    static_cast<std::uint64_t>(word_index),
                                    &value),
              "gbits word query");
        SEXP out = PROTECT(Rf_allocVector(RAWSXP, 8));
        for (unsigned int byte = 0; byte < 8u; ++byte) {
            RAW(out)[byte] = static_cast<Rbyte>((value >> (byte * 8u)) & 0xffu);
        }
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("gbits word query: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_copy_interval(
    SEXP destination_pointer, SEXP destination_individual_sexp,
    SEXP source_pointer, SEXP source_individual_sexp,
    SEXP first_marker_sexp, SEXP last_marker_sexp) {
    try {
        Packed* destination = require_packed(destination_pointer);
        Packed* source = require_packed(source_pointer);
        require_same_backend(destination, source);
        check(destination->backend,
              destination->backend->copy_interval(
                  destination->handle,
                  static_cast<std::uint64_t>(scalar_int(
                      destination_individual_sexp, "destination individual")),
                  source->handle,
                  static_cast<std::uint64_t>(scalar_int(
                      source_individual_sexp, "source individual")),
                  static_cast<std::uint64_t>(scalar_int(first_marker_sexp,
                                                        "first marker")),
                  static_cast<std::uint64_t>(scalar_int(last_marker_sexp,
                                                        "last marker"))),
              "gbits interval copy");
        return destination_pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits interval copy: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_copy_filtered(
    SEXP destination_pointer, SEXP destination_individual_sexp,
    SEXP source_pointer, SEXP source_individual_sexp,
    SEXP first_marker_sexp, SEXP last_marker_sexp, SEXP age_sexp,
    SEXP mutation_age) {
    try {
        Packed* destination = require_packed(destination_pointer);
        Packed* source = require_packed(source_pointer);
        require_same_backend(destination, source);
        if (TYPEOF(age_sexp) != REALSXP || XLENGTH(age_sexp) != 1 ||
            !R_FINITE(REAL(age_sexp)[0]) || TYPEOF(mutation_age) != REALSXP) {
            fail("filtered copy requires numeric age inputs");
        }
        check(destination->backend,
              destination->backend->copy_filtered(
                  destination->handle,
                  static_cast<std::uint64_t>(scalar_int(
                      destination_individual_sexp, "destination individual")),
                  source->handle,
                  static_cast<std::uint64_t>(scalar_int(
                      source_individual_sexp, "source individual")),
                  static_cast<std::uint64_t>(scalar_int(first_marker_sexp,
                                                        "first marker")),
                  static_cast<std::uint64_t>(scalar_int(last_marker_sexp,
                                                        "last marker")),
                  REAL(age_sexp)[0], REAL(mutation_age),
                  static_cast<std::uint64_t>(XLENGTH(mutation_age))),
              "gbits filtered copy");
        return destination_pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits filtered copy: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_make_gamete(
    SEXP destination_pointer, SEXP destination_individual_sexp,
    SEXP parent_h1_pointer, SEXP parent_h2_pointer,
    SEXP parent_individual_sexp, SEXP starting_haplotype_sexp,
    SEXP boundaries) {
    try {
        Packed* destination = require_packed(destination_pointer);
        Packed* h1 = require_packed(parent_h1_pointer);
        Packed* h2 = require_packed(parent_h2_pointer);
        require_same_backend(destination, h1);
        require_same_backend(destination, h2);
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
        check(destination->backend,
              destination->backend->make_gamete(
                  destination->handle,
                  static_cast<std::uint64_t>(scalar_int(
                      destination_individual_sexp, "destination individual")),
                  h1->handle, h2->handle,
                  static_cast<std::uint64_t>(scalar_int(
                      parent_individual_sexp, "parent individual")),
                  static_cast<std::uint32_t>(starting - 1),
                  converted.empty() ? nullptr : converted.data(),
                  static_cast<std::uint64_t>(converted.size())),
              "gbits gamete construction");
        return destination_pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits gamete construction: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_decode_genotypes(SEXP h1_pointer,
                                               SEXP h2_pointer) {
    try {
        Packed* h1 = require_packed(h1_pointer);
        Packed* h2 = require_packed(h2_pointer);
        require_same_backend(h1, h2);
        if (h1->individuals != h2->individuals || h1->markers != h2->markers ||
            h1->individuals > static_cast<std::uint64_t>(
                                  std::numeric_limits<int>::max()) ||
            h1->markers > static_cast<std::uint64_t>(
                              std::numeric_limits<int>::max())) {
            fail("packed phases have incompatible R dimensions");
        }
        SEXP out = PROTECT(Rf_allocMatrix(
            RAWSXP, static_cast<int>(h1->individuals),
            static_cast<int>(h1->markers)));
        check(h1->backend,
              h1->backend->decode_genotypes(
                  h1->handle, h2->handle, RAW(out),
                  static_cast<std::uint64_t>(XLENGTH(out)), h1->individuals),
              "gbits genotype decoding");
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("gbits genotype decoding: %s", ex.what());
    }
    return R_NilValue;
}
