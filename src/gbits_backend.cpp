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
struct BedSinkInfo;
struct HapSinkInfo;
struct HapChromosomeInfo;

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
    status_t (*bed_open)(const char*, std::uint64_t, std::uint64_t, handle_t**);
    status_t (*bed_close)(handle_t*);
    status_t (*bed_read_variant)(handle_t*, std::uint64_t, std::int8_t*,
                                 std::uint64_t);
    status_t (*bed_sink_create)(const char*, std::uint64_t, std::uint32_t,
                                std::uint64_t, handle_t**);
    status_t (*bed_sink_append)(handle_t*, const handle_t*, const handle_t*);
    status_t (*bed_sink_finalize)(handle_t*);
    status_t (*bed_sink_info)(const handle_t*, BedSinkInfo*);
    status_t (*bed_sink_close)(handle_t*);
    status_t (*hap_sink_create)(const char*, std::uint64_t, std::uint32_t,
                                handle_t**);
    status_t (*hap_sink_append)(handle_t*, const handle_t*, const handle_t*);
    status_t (*hap_sink_finalize)(handle_t*);
    status_t (*hap_sink_info)(const handle_t*, HapSinkInfo*);
    status_t (*hap_sink_close)(handle_t*);
    status_t (*hap_reader_open)(const char*, handle_t**);
    status_t (*hap_reader_close)(handle_t*);
    status_t (*hap_reader_dimensions)(const handle_t*, std::uint64_t*,
                                      std::uint64_t*, std::uint64_t*);
    status_t (*hap_reader_chromosome_info)(const handle_t*, std::uint64_t,
                                           HapChromosomeInfo*);
    status_t (*hap_reader_load)(const handle_t*, std::uint64_t, handle_t**,
                                handle_t**);
    std::string version;
};

struct Packed {
    Backend* backend;
    handle_t* handle;
    std::uint64_t individuals;
    std::uint64_t markers;
};

struct BedSinkInfo {
    std::uint64_t individual_count;
    std::uint64_t variant_count;
    std::uint64_t bytes_written;
    std::uint64_t conversion_buffer_bytes;
    std::uint64_t lifecycle_object_bytes;
    int state;
};

struct BedSink {
    Backend* backend;
    handle_t* handle;
};

struct HapSinkInfo {
    std::uint64_t individual_count;
    std::uint64_t marker_count;
    std::uint64_t chromosome_count;
    std::uint64_t bytes_written;
    int state;
};

struct HapChromosomeInfo {
    std::uint64_t global_start_marker;
    std::uint64_t marker_count;
    std::uint64_t h1_offset;
    std::uint64_t h2_offset;
    std::uint64_t bytes_per_phase;
};

struct HapSink { Backend* backend; handle_t* handle; };
struct HapReader {
    Backend* backend;
    handle_t* handle;
    std::uint64_t individuals;
    std::uint64_t markers;
    std::uint64_t chromosomes;
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

void bed_sink_finalizer(SEXP pointer) {
    BedSink* sink = static_cast<BedSink*>(R_ExternalPtrAddr(pointer));
    if (sink != nullptr) {
        if (sink->handle != nullptr && sink->backend != nullptr) {
            (void)sink->backend->bed_sink_close(sink->handle);
        }
        delete sink;
    }
    R_ClearExternalPtr(pointer);
}

void hap_sink_finalizer(SEXP pointer) {
    HapSink* sink = static_cast<HapSink*>(R_ExternalPtrAddr(pointer));
    if (sink != nullptr) {
        if (sink->handle != nullptr) (void)sink->backend->hap_sink_close(sink->handle);
        delete sink;
    }
    R_ClearExternalPtr(pointer);
}

void hap_reader_finalizer(SEXP pointer) {
    HapReader* reader = static_cast<HapReader*>(R_ExternalPtrAddr(pointer));
    if (reader != nullptr) {
        if (reader->handle != nullptr) (void)reader->backend->hap_reader_close(reader->handle);
        delete reader;
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

BedSink* require_bed_sink(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP) fail("BED sink is invalid");
    BedSink* sink = static_cast<BedSink*>(R_ExternalPtrAddr(pointer));
    if (sink == nullptr || sink->handle == nullptr) {
        fail("BED sink has been released");
    }
    return sink;
}

HapSink* require_hap_sink(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP) fail("HAP sink is invalid");
    HapSink* sink = static_cast<HapSink*>(R_ExternalPtrAddr(pointer));
    if (sink == nullptr || sink->handle == nullptr) fail("HAP sink has been released");
    return sink;
}

HapReader* require_hap_reader(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP) fail("HAP reader is invalid");
    HapReader* reader = static_cast<HapReader*>(R_ExternalPtrAddr(pointer));
    if (reader == nullptr || reader->handle == nullptr) fail("HAP reader has been released");
    return reader;
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

double exact_r_number(std::uint64_t value, const char* field) {
    if (value > 9007199254740991ULL) {
        fail(std::string(field) + " exceeds exact R numeric representation");
    }
    return static_cast<double>(value);
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
            backend->bed_open = symbol<decltype(backend->bed_open)>(
                symbols, "gbits_bed_open");
            backend->bed_close = symbol<decltype(backend->bed_close)>(
                symbols, "gbits_bed_close");
            backend->bed_read_variant = symbol<decltype(backend->bed_read_variant)>(
                symbols, "gbits_bed_read_variant");
            backend->bed_sink_create = symbol<decltype(backend->bed_sink_create)>(
                symbols, "gbits_bed_sink_create");
            backend->bed_sink_append = symbol<decltype(backend->bed_sink_append)>(
                symbols, "gbits_bed_sink_append_phased");
            backend->bed_sink_finalize = symbol<decltype(backend->bed_sink_finalize)>(
                symbols, "gbits_bed_sink_finalize");
            backend->bed_sink_info = symbol<decltype(backend->bed_sink_info)>(
                symbols, "gbits_bed_sink_get_info");
            backend->bed_sink_close = symbol<decltype(backend->bed_sink_close)>(
                symbols, "gbits_bed_sink_close");
            backend->hap_sink_create = symbol<decltype(backend->hap_sink_create)>(
                symbols, "gbits_hap_sink_create");
            backend->hap_sink_append = symbol<decltype(backend->hap_sink_append)>(
                symbols, "gbits_hap_sink_append_phased");
            backend->hap_sink_finalize = symbol<decltype(backend->hap_sink_finalize)>(
                symbols, "gbits_hap_sink_finalize");
            backend->hap_sink_info = symbol<decltype(backend->hap_sink_info)>(
                symbols, "gbits_hap_sink_get_info");
            backend->hap_sink_close = symbol<decltype(backend->hap_sink_close)>(
                symbols, "gbits_hap_sink_close");
            backend->hap_reader_open = symbol<decltype(backend->hap_reader_open)>(
                symbols, "gbits_hap_reader_open");
            backend->hap_reader_close = symbol<decltype(backend->hap_reader_close)>(
                symbols, "gbits_hap_reader_close");
            backend->hap_reader_dimensions = symbol<decltype(backend->hap_reader_dimensions)>(
                symbols, "gbits_hap_reader_dimensions");
            backend->hap_reader_chromosome_info = symbol<decltype(backend->hap_reader_chromosome_info)>(
                symbols, "gbits_hap_reader_chromosome_info");
            backend->hap_reader_load = symbol<decltype(backend->hap_reader_load)>(
                symbols, "gbits_hap_reader_load_chromosome");
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

extern "C" SEXP C_gsim_gbits_close(SEXP pointer) {
    try {
        Packed* packed = require_packed(pointer);
        Backend* backend = packed->backend;
        handle_t* handle = packed->handle;
        packed->handle = nullptr;
        delete packed;
        R_ClearExternalPtr(pointer);
        check(backend, backend->close(handle), "gbits packed haplotype close");
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("gbits packed haplotype close: %s", ex.what());
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

extern "C" SEXP C_gsim_gbits_copy_filtered_counts(
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
        const int destination_individual = scalar_int(
            destination_individual_sexp, "destination individual");
        const int source_individual = scalar_int(
            source_individual_sexp, "source individual");
        const int first_marker = scalar_int(first_marker_sexp, "first marker");
        const int last_marker = scalar_int(last_marker_sexp, "last marker");
        check(destination->backend,
              destination->backend->copy_filtered(
                  destination->handle,
                  static_cast<std::uint64_t>(destination_individual),
                  source->handle, static_cast<std::uint64_t>(source_individual),
                  static_cast<std::uint64_t>(first_marker),
                  static_cast<std::uint64_t>(last_marker), REAL(age_sexp)[0],
                  REAL(mutation_age),
                  static_cast<std::uint64_t>(XLENGTH(mutation_age))),
              "gbits filtered copy");
        int copied = 0;
        int retained = 0;
        for (int marker = first_marker;; ++marker) {
            std::uint8_t allele = 0u;
            check(source->backend,
                  source->backend->allele(
                      source->handle,
                      static_cast<std::uint64_t>(source_individual),
                      static_cast<std::uint64_t>(marker), &allele),
                  "gbits packed allele query");
            copied += allele == 1u ? 1 : 0;
            retained += allele == 1u &&
                        REAL(age_sexp)[0] < REAL(mutation_age)[marker]
                            ? 1 : 0;
            if (marker == last_marker) break;
        }
        SEXP out = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(out)[0] = copied;
        INTEGER(out)[1] = retained;
        set_names(out, {"copied_alternative", "retained_alternative"});
        UNPROTECT(1);
        return out;
    } catch (const std::exception& ex) {
        Rf_error("gbits filtered copy with audit: %s", ex.what());
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

extern "C" SEXP C_gsim_gbits_bed_sink_create(
    SEXP backend_pointer, SEXP path_sexp, SEXP individuals_sexp,
    SEXP overwrite_sexp, SEXP buffer_variants_sexp) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const std::string path = scalar_utf8(path_sexp, "BED path");
        const int individuals = scalar_int(individuals_sexp, "individuals", 1);
        const int buffer_variants =
            scalar_int(buffer_variants_sexp, "buffer_variants", 1);
        const bool overwrite = scalar_bool(overwrite_sexp, "overwrite");
        handle_t* handle = nullptr;
        check(backend,
              backend->bed_sink_create(
                  path.c_str(), static_cast<std::uint64_t>(individuals),
                  overwrite ? 1u : 0u,
                  static_cast<std::uint64_t>(buffer_variants), &handle),
              "gbits BED sink creation");
        BedSink* sink = nullptr;
        try {
            sink = new BedSink{backend, handle};
        } catch (...) {
            (void)backend->bed_sink_close(handle);
            throw;
        }
        SEXP pointer = PROTECT(R_MakeExternalPtr(sink, R_NilValue,
                                                 backend_pointer));
        R_RegisterCFinalizerEx(pointer, bed_sink_finalizer, TRUE);
        UNPROTECT(1);
        return pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits BED sink creation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_bed_sink_append(
    SEXP sink_pointer, SEXP h1_pointer, SEXP h2_pointer) {
    try {
        BedSink* sink = require_bed_sink(sink_pointer);
        Packed* h1 = require_packed(h1_pointer);
        Packed* h2 = require_packed(h2_pointer);
        if (sink->backend != h1->backend || sink->backend != h2->backend) {
            fail("BED sink and packed phases originate from different gbits backends");
        }
        check(sink->backend,
              sink->backend->bed_sink_append(sink->handle, h1->handle,
                                             h2->handle),
              "gbits BED chromosome append");
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits BED chromosome append: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_bed_sink_finalize(SEXP sink_pointer) {
    try {
        BedSink* sink = require_bed_sink(sink_pointer);
        check(sink->backend, sink->backend->bed_sink_finalize(sink->handle),
              "gbits BED sink finalization");
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits BED sink finalization: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_bed_sink_cancel(SEXP sink_pointer) {
    try {
        BedSink* sink = require_bed_sink(sink_pointer);
        check(sink->backend, sink->backend->bed_sink_close(sink->handle),
              "gbits BED sink cancellation");
        sink->handle = nullptr;
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("gbits BED sink cancellation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_bed_sink_info(SEXP sink_pointer) {
    try {
        BedSink* sink = require_bed_sink(sink_pointer);
        BedSinkInfo info{};
        check(sink->backend,
              sink->backend->bed_sink_info(sink->handle, &info),
              "gbits BED sink information");
        SEXP values = PROTECT(Rf_allocVector(REALSXP, 5));
        REAL(values)[0] = static_cast<double>(info.individual_count);
        REAL(values)[1] = static_cast<double>(info.variant_count);
        REAL(values)[2] = static_cast<double>(info.bytes_written);
        REAL(values)[3] = static_cast<double>(info.conversion_buffer_bytes);
        REAL(values)[4] = static_cast<double>(info.lifecycle_object_bytes);
        set_names(values, {"individual_count", "variant_count",
                           "bytes_written", "conversion_buffer_bytes",
                           "lifecycle_object_bytes"});
        const char* state = info.state == 0 ? "open"
                            : info.state == 1 ? "finalized"
                                              : "failed";
        SEXP state_value = PROTECT(Rf_mkString(state));
        Rf_setAttrib(values, Rf_install("state"), state_value);
        UNPROTECT(2);
        return values;
    } catch (const std::exception& ex) {
        Rf_error("gbits BED sink information: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_hap_sink_create(
    SEXP backend_pointer, SEXP path_sexp, SEXP individuals_sexp,
    SEXP overwrite_sexp) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const std::string path = scalar_utf8(path_sexp, "HAP path");
        const int individuals = scalar_int(individuals_sexp, "individuals", 1);
        const bool overwrite = scalar_bool(overwrite_sexp, "overwrite");
        handle_t* handle = nullptr;
        check(backend, backend->hap_sink_create(
              path.c_str(), static_cast<std::uint64_t>(individuals),
              overwrite ? 1u : 0u, &handle), "gbits HAP sink creation");
        HapSink* sink = nullptr;
        try { sink = new HapSink{backend, handle}; }
        catch (...) { (void)backend->hap_sink_close(handle); throw; }
        SEXP pointer = PROTECT(R_MakeExternalPtr(sink, R_NilValue, backend_pointer));
        R_RegisterCFinalizerEx(pointer, hap_sink_finalizer, TRUE);
        UNPROTECT(1);
        return pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits HAP sink creation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_hap_sink_append(
    SEXP sink_pointer, SEXP h1_pointer, SEXP h2_pointer) {
    try {
        HapSink* sink = require_hap_sink(sink_pointer);
        Packed* h1 = require_packed(h1_pointer);
        Packed* h2 = require_packed(h2_pointer);
        if (sink->backend != h1->backend || sink->backend != h2->backend) {
            fail("HAP sink and phases originate from different gbits backends");
        }
        check(sink->backend, sink->backend->hap_sink_append(
              sink->handle, h1->handle, h2->handle), "gbits HAP append");
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits HAP append: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_hap_sink_finalize(SEXP sink_pointer) {
    try {
        HapSink* sink = require_hap_sink(sink_pointer);
        check(sink->backend, sink->backend->hap_sink_finalize(sink->handle),
              "gbits HAP finalization");
        return sink_pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits HAP finalization: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_hap_sink_cancel(SEXP sink_pointer) {
    try {
        HapSink* sink = require_hap_sink(sink_pointer);
        check(sink->backend, sink->backend->hap_sink_close(sink->handle),
              "gbits HAP cancellation");
        sink->handle = nullptr;
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("gbits HAP cancellation: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_hap_sink_info(SEXP sink_pointer) {
    try {
        HapSink* sink = require_hap_sink(sink_pointer);
        HapSinkInfo info{};
        check(sink->backend, sink->backend->hap_sink_info(sink->handle, &info),
              "gbits HAP sink information");
        SEXP values = PROTECT(Rf_allocVector(REALSXP, 4));
        REAL(values)[0] = exact_r_number(info.individual_count, "HAP individual count");
        REAL(values)[1] = exact_r_number(info.marker_count, "HAP marker count");
        REAL(values)[2] = exact_r_number(info.chromosome_count, "HAP chromosome count");
        REAL(values)[3] = exact_r_number(info.bytes_written, "HAP byte count");
        set_names(values, {"individual_count", "marker_count",
                           "chromosome_count", "bytes_written"});
        const char* state = info.state == 0 ? "open" :
                            info.state == 1 ? "finalized" : "failed";
        SEXP state_value = PROTECT(Rf_mkString(state));
        Rf_setAttrib(values, Rf_install("state"), state_value);
        UNPROTECT(2);
        return values;
    } catch (const std::exception& ex) {
        Rf_error("gbits HAP sink information: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_hap_reader_open(
    SEXP backend_pointer, SEXP path_sexp) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const std::string path = scalar_utf8(path_sexp, "HAP path");
        handle_t* handle = nullptr;
        check(backend, backend->hap_reader_open(path.c_str(), &handle),
              "gbits HAP reader open");
        std::uint64_t individuals = 0u, markers = 0u, chromosomes = 0u;
        try {
            check(backend, backend->hap_reader_dimensions(
                  handle, &individuals, &markers, &chromosomes),
                  "gbits HAP reader dimensions");
        } catch (...) { (void)backend->hap_reader_close(handle); throw; }
        HapReader* reader = nullptr;
        try { reader = new HapReader{backend, handle, individuals, markers, chromosomes}; }
        catch (...) { (void)backend->hap_reader_close(handle); throw; }
        SEXP pointer = PROTECT(R_MakeExternalPtr(reader, R_NilValue, backend_pointer));
        R_RegisterCFinalizerEx(pointer, hap_reader_finalizer, TRUE);
        UNPROTECT(1);
        return pointer;
    } catch (const std::exception& ex) {
        Rf_error("gbits HAP reader open: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_hap_reader_close(SEXP reader_pointer) {
    try {
        HapReader* reader = require_hap_reader(reader_pointer);
        check(reader->backend, reader->backend->hap_reader_close(reader->handle),
              "gbits HAP reader close");
        reader->handle = nullptr;
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("gbits HAP reader close: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_hap_reader_info(SEXP reader_pointer) {
    try {
        HapReader* reader = require_hap_reader(reader_pointer);
        if (reader->chromosomes > static_cast<std::uint64_t>(
                std::numeric_limits<int>::max())) fail("too many HAP chromosomes for R");
        const int count = static_cast<int>(reader->chromosomes);
        SEXP ranges = PROTECT(Rf_allocMatrix(REALSXP, count, 5));
        for (int i = 0; i < count; ++i) {
            HapChromosomeInfo info{};
            check(reader->backend, reader->backend->hap_reader_chromosome_info(
                  reader->handle, static_cast<std::uint64_t>(i), &info),
                  "gbits HAP chromosome information");
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
        SET_VECTOR_ELT(result, 0, Rf_ScalarReal(exact_r_number(reader->individuals, "HAP individual count")));
        SET_VECTOR_ELT(result, 1, Rf_ScalarReal(exact_r_number(reader->markers, "HAP marker count")));
        SET_VECTOR_ELT(result, 2, Rf_ScalarReal(exact_r_number(reader->chromosomes, "HAP chromosome count")));
        SET_VECTOR_ELT(result, 3, ranges);
        set_names(result, {"individual_count", "marker_count", "chromosome_count", "ranges"});
        UNPROTECT(4);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("gbits HAP reader information: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_hap_reader_load(
    SEXP reader_pointer, SEXP chromosome_sexp) {
    try {
        HapReader* reader = require_hap_reader(reader_pointer);
        const int chromosome = scalar_int(chromosome_sexp, "chromosome", 1);
        if (static_cast<std::uint64_t>(chromosome) > reader->chromosomes) {
            fail("HAP chromosome is out of range");
        }
        handle_t* h1 = nullptr;
        handle_t* h2 = nullptr;
        check(reader->backend, reader->backend->hap_reader_load(
              reader->handle, static_cast<std::uint64_t>(chromosome - 1), &h1, &h2),
              "gbits HAP chromosome load");
        HapChromosomeInfo info{};
        check(reader->backend, reader->backend->hap_reader_chromosome_info(
              reader->handle, static_cast<std::uint64_t>(chromosome - 1), &info),
              "gbits HAP loaded chromosome information");
        SEXP backend_pointer = R_ExternalPtrProtected(reader_pointer);
        SEXP first = PROTECT(make_packed(backend_pointer, h1, reader->individuals,
                                         info.marker_count));
        SEXP second = PROTECT(make_packed(backend_pointer, h2, reader->individuals,
                                          info.marker_count));
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(result, 0, first);
        SET_VECTOR_ELT(result, 1, second);
        set_names(result, {"h1", "h2"});
        UNPROTECT(3);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("gbits HAP chromosome load: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_gbits_bed_read_all(
    SEXP backend_pointer, SEXP path_sexp, SEXP individuals_sexp,
    SEXP variants_sexp) {
    try {
        Backend* backend = require_backend(backend_pointer);
        const std::string path = scalar_utf8(path_sexp, "BED path");
        const int individuals = scalar_int(individuals_sexp, "individuals", 1);
        const int variants = scalar_int(variants_sexp, "variants", 1);
        handle_t* reader = nullptr;
        check(backend,
              backend->bed_open(path.c_str(),
                                static_cast<std::uint64_t>(individuals),
                                static_cast<std::uint64_t>(variants), &reader),
              "gbits BED reader open");
        SEXP output = PROTECT(Rf_allocMatrix(INTSXP, individuals, variants));
        std::vector<std::int8_t> record(static_cast<std::size_t>(individuals));
        try {
            for (int marker = 0; marker < variants; ++marker) {
                check(backend,
                      backend->bed_read_variant(
                          reader, static_cast<std::uint64_t>(marker),
                          record.data(), static_cast<std::uint64_t>(individuals)),
                      "gbits BED reader decode");
                for (int individual = 0; individual < individuals; ++individual) {
                    INTEGER(output)[individual + individuals * marker] =
                        static_cast<int>(record[static_cast<std::size_t>(individual)]);
                }
            }
        } catch (...) {
            (void)backend->bed_close(reader);
            UNPROTECT(1);
            throw;
        }
        check(backend, backend->bed_close(reader), "gbits BED reader close");
        UNPROTECT(1);
        return output;
    } catch (const std::exception& ex) {
        Rf_error("gbits BED validation decode: %s", ex.what());
    }
    return R_NilValue;
}
