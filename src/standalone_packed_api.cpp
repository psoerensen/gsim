#include "standalone_packed_api.h"

#include "standalone_bed.h"
#include "standalone_error.h"
#include "standalone_hap.h"
#include "standalone_packed.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace gsim::native::api {
namespace {

struct PackedHandle {
    explicit PackedHandle(PhasedHaplotypeMatrix value) : value(std::move(value)) {}
    PhasedHaplotypeMatrix value;
};
struct BedSinkHandle {
    BedSinkHandle(std::string path, std::uint64_t n, bool overwrite,
                  std::uint64_t capacity)
        : value(std::move(path), n, overwrite, capacity) {}
    PhasedBedSink value;
};
struct HapSinkHandle {
    HapSinkHandle(std::string path, std::uint64_t n, bool overwrite)
        : value(std::move(path), n, overwrite) {}
    PhasedHapSink value;
};
struct HapReaderHandle {
    explicit HapReaderHandle(std::string path) : value(std::move(path)) {}
    PhasedHapReader value;
};
struct BedReaderHandle {
    BedReaderHandle(const char* path, std::uint64_t samples,
                    std::uint64_t variants)
        : input(std::filesystem::u8path(path), std::ios::binary), n(samples),
          m(variants), bytes_per_variant(samples / 4u + (samples % 4u != 0u)) {
        if (n == 0u || m == 0u) invalid_argument("BED dimensions must be positive");
        if (!input.is_open()) throw Error(status_io_error, "cannot open BED file");
        std::uint8_t header[3]{};
        input.read(reinterpret_cast<char*>(header), 3);
        if (!input || header[0] != 0x6cu || header[1] != 0x1bu ||
            header[2] != 0x01u) {
            throw Error(status_bed_format_error, "invalid SNP-major BED header");
        }
        if (m > (std::numeric_limits<std::uint64_t>::max() - 3u) /
                    bytes_per_variant) {
            invalid_argument("BED size overflows uint64_t");
        }
        const std::uint64_t expected = 3u + m * bytes_per_variant;
        input.seekg(0, std::ios::end);
        const auto size = input.tellg();
        if (size < 0 || static_cast<std::uint64_t>(size) != expected) {
            throw Error(status_bed_format_error, "BED file size does not match dimensions");
        }
    }
    std::ifstream input;
    std::uint64_t n, m, bytes_per_variant;
};

thread_local std::string error_message;

template<class Function> status_t protect(Function&& function) noexcept {
    try {
        function();
        error_message.clear();
        return status_success;
    } catch (const Error& error) {
        error_message = error.what();
        return static_cast<int>(error.status());
    } catch (const std::bad_alloc&) {
        error_message = "memory allocation failed";
        return status_internal_error;
    } catch (const std::exception& error) {
        error_message = error.what();
        return status_internal_error;
    } catch (...) {
        error_message = "unknown internal error";
        return status_internal_error;
    }
}

template<class T> T& required(handle_t* value, const char* name) {
    if (value == nullptr) invalid_argument(std::string(name) + " handle is null");
    return *static_cast<T*>(value);
}
template<class T> const T& required(const handle_t* value, const char* name) {
    if (value == nullptr) invalid_argument(std::string(name) + " handle is null");
    return *static_cast<const T*>(value);
}
const PhasedHaplotypeMatrix& packed(const handle_t* value) {
    return required<PackedHandle>(value, "packed haplotype").value;
}
PhasedHaplotypeMatrix& packed(handle_t* value) {
    return required<PackedHandle>(value, "packed haplotype").value;
}

} // namespace

std::uint32_t abi_version() noexcept { return 1u; }
const char* library_version() noexcept { return "gsim-native-1"; }
const char* last_error() noexcept { return error_message.c_str(); }

status_t create_zero(std::uint64_t n, std::uint64_t m, handle_t** out) {
    if (out) *out = nullptr;
    return protect([&] { if (!out) invalid_argument("output handle is null");
        *out = new PackedHandle(PhasedHaplotypeMatrix(n, m)); });
}
status_t create_values(std::uint64_t n, std::uint64_t m,
                       const std::uint8_t* values, std::uint64_t length,
                       std::uint64_t ld, handle_t** out) {
    if (out) *out = nullptr;
    return protect([&] { if (!out) invalid_argument("output handle is null");
        *out = new PackedHandle(PhasedHaplotypeMatrix::from_values(
            n, m, values, length, ld)); });
}
status_t close(handle_t* value) { return protect([&] { delete &required<PackedHandle>(value, "packed haplotype"); }); }
status_t individual_count(const handle_t* value, std::uint64_t* out) {
    return protect([&] { if (!out) invalid_argument("count output is null"); *out = packed(value).individual_count(); });
}
status_t marker_count(const handle_t* value, std::uint64_t* out) {
    return protect([&] { if (!out) invalid_argument("count output is null"); *out = packed(value).marker_count(); });
}
status_t words_per_marker(const handle_t* value, std::uint64_t* out) {
    return protect([&] { if (!out) invalid_argument("count output is null"); *out = packed(value).words_per_marker(); });
}
status_t storage_bytes(const handle_t* value, std::uint64_t* out) {
    return protect([&] { if (!out) invalid_argument("byte output is null"); *out = packed(value).storage_bytes(); });
}
status_t word(const handle_t* value, std::uint64_t marker, std::uint64_t index,
              std::uint64_t* out) {
    return protect([&] { if (!out) invalid_argument("word output is null"); *out = packed(value).word(marker, index); });
}
status_t allele(const handle_t* value, std::uint64_t individual,
                std::uint64_t marker, std::uint8_t* out) {
    return protect([&] { if (!out) invalid_argument("allele output is null"); *out = packed(value).allele(individual, marker); });
}
status_t set_allele(handle_t* value, std::uint64_t individual,
                    std::uint64_t marker, std::uint8_t x) {
    return protect([&] { packed(value).set_allele(individual, marker, x); });
}
status_t unpack(const handle_t* value, std::uint8_t* out, std::uint64_t length,
                std::uint64_t ld) {
    return protect([&] { packed(value).unpack(out, length, ld); });
}
status_t copy_interval(handle_t* destination, std::uint64_t di,
                       const handle_t* source, std::uint64_t si,
                       std::uint64_t first, std::uint64_t last) {
    return protect([&] { packed(destination).copy_interval(di, packed(source), si, first, last); });
}
status_t copy_filtered(handle_t* destination, std::uint64_t di,
                       const handle_t* source, std::uint64_t si,
                       std::uint64_t first, std::uint64_t last, double age,
                       const double* mutation, std::uint64_t count) {
    return protect([&] { packed(destination).copy_filtered_segment(
        di, packed(source), si, first, last, age, mutation, count); });
}
status_t materialize_founders(
    handle_t* destination_h1, handle_t* destination_h2,
    const handle_t* reference_h1, const handle_t* reference_h2,
    const std::uint64_t* destination, const std::uint32_t* phase,
    const std::uint64_t* donor, const std::uint64_t* first,
    const std::uint64_t* last, const double* age, std::uint64_t event_count,
    const double* mutation, std::uint64_t mutation_count,
    std::uint32_t requested_threads, std::uint64_t* copied,
    std::uint64_t* retained) {
    return protect([&] {
        auto& out_h1 = packed(destination_h1);
        auto& out_h2 = packed(destination_h2);
        const auto& ref_h1 = packed(reference_h1);
        const auto& ref_h2 = packed(reference_h2);
        if (out_h1.individual_count() != out_h2.individual_count() ||
            out_h1.marker_count() != out_h2.marker_count() ||
            ref_h1.individual_count() != ref_h2.individual_count() ||
            ref_h1.marker_count() != ref_h2.marker_count() ||
            out_h1.marker_count() != ref_h1.marker_count()) {
            invalid_argument("founder materialization requires compatible phases");
        }
        if (requested_threads == 0u || mutation == nullptr ||
            mutation_count != out_h1.marker_count() ||
            (event_count != 0u &&
             (destination == nullptr || phase == nullptr || donor == nullptr ||
              first == nullptr || last == nullptr || age == nullptr))) {
            invalid_argument("founder materialization arguments are invalid");
        }
        for (std::uint64_t marker = 0u; marker < mutation_count; ++marker) {
            if (!std::isfinite(mutation[static_cast<std::size_t>(marker)]) ||
                mutation[static_cast<std::size_t>(marker)] < 0.0) {
                invalid_argument("mutation ages must be finite and nonnegative");
            }
        }
        std::uint64_t previous = 0u;
        for (std::uint64_t i = 0u; i < event_count; ++i) {
            if (destination[i] >= out_h1.individual_count() || phase[i] > 1u ||
                donor[i] >= ref_h1.individual_count() || first[i] > last[i] ||
                last[i] >= out_h1.marker_count() || !std::isfinite(age[i]) ||
                age[i] < 0.0 || (i != 0u && destination[i] < previous)) {
                invalid_argument("founder event plan is invalid or unsorted");
            }
            previous = destination[i];
        }
        if (event_count == 0u) return;
        const std::uint64_t words = out_h1.words_per_marker();
        const std::uint32_t worker_count = static_cast<std::uint32_t>(
            std::min<std::uint64_t>(requested_threads, words));
        std::vector<std::exception_ptr> failures(worker_count);
        auto work = [&](std::uint32_t worker) {
            try {
                const std::uint64_t first_word = words * worker / worker_count;
                const std::uint64_t last_word = words * (worker + 1u) / worker_count;
                const std::uint64_t first_individual = first_word * 64u;
                const std::uint64_t last_individual = std::min(
                    out_h1.individual_count(), last_word * 64u);
                const std::uint64_t* begin = std::lower_bound(
                    destination, destination + event_count, first_individual);
                const std::uint64_t* end = std::lower_bound(
                    begin, destination + event_count, last_individual);
                for (const std::uint64_t* item = begin; item != end; ++item) {
                    const std::uint64_t i = static_cast<std::uint64_t>(item - destination);
                    auto& output = phase[i] == 0u ? out_h1 : out_h2;
                    const auto& source = phase[i] == 0u ? ref_h1 : ref_h2;
                    const auto counts = output.copy_filtered_segment_counts(
                        destination[i], source, donor[i], first[i], last[i],
                        age[i], mutation, mutation_count);
                    if (copied != nullptr) copied[i] = counts.first;
                    if (retained != nullptr) retained[i] = counts.second;
                }
            } catch (...) {
                failures[worker] = std::current_exception();
            }
        };
        std::vector<std::thread> workers;
        workers.reserve(worker_count > 0u ? worker_count - 1u : 0u);
        try {
            for (std::uint32_t worker = 1u; worker < worker_count; ++worker) {
                workers.emplace_back(work, worker);
            }
            work(0u);
        } catch (...) {
            for (auto& worker : workers) if (worker.joinable()) worker.join();
            throw;
        }
        for (auto& worker : workers) worker.join();
        for (const auto& failure : failures) if (failure) std::rethrow_exception(failure);
    });
}
status_t make_gamete(handle_t* destination, std::uint64_t di,
                     const handle_t* h1, const handle_t* h2, std::uint64_t pi,
                     std::uint32_t starting, const std::uint64_t* boundaries,
                     std::uint64_t count) {
    return protect([&] { packed(destination).make_gamete(
        di, packed(h1), packed(h2), pi, starting, boundaries, count); });
}
status_t decode_genotypes(const handle_t* h1, const handle_t* h2,
                          std::uint8_t* out, std::uint64_t length,
                          std::uint64_t ld) {
    return protect([&] { packed(h1).decode_genotypes(packed(h2), out, length, ld); });
}

status_t bed_open(const char* path, std::uint64_t n, std::uint64_t m,
                  handle_t** out) {
    if (out) *out = nullptr;
    return protect([&] { if (!out || !path || !*path) invalid_argument("BED path/output is empty"); *out = new BedReaderHandle(path, n, m); });
}
status_t bed_close(handle_t* value) { return protect([&] { delete &required<BedReaderHandle>(value, "BED reader"); }); }
status_t bed_read_variant(handle_t* value, std::uint64_t marker,
                          std::int8_t* out, std::uint64_t length) {
    return protect([&] {
        auto& reader = required<BedReaderHandle>(value, "BED reader");
        if (!out || length != reader.n || marker >= reader.m) invalid_argument("BED variant request is invalid");
        std::vector<std::uint8_t> bytes(static_cast<std::size_t>(reader.bytes_per_variant));
        const std::uint64_t offset = 3u + marker * reader.bytes_per_variant;
        reader.input.clear(); reader.input.seekg(static_cast<std::streamoff>(offset));
        reader.input.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        if (!reader.input) throw Error(status_io_error, "cannot read BED variant");
        for (std::uint64_t i = 0; i < reader.n; ++i) {
            const std::uint8_t code = (bytes[static_cast<std::size_t>(i / 4u)] >> (2u * (i % 4u))) & 3u;
            static const std::int8_t dosage[4] = {2, -1, 1, 0};
            out[static_cast<std::size_t>(i)] = dosage[code];
        }
    });
}

status_t bed_sink_create(const char* path, std::uint64_t n,
                         std::uint32_t overwrite, std::uint64_t capacity,
                         handle_t** out) {
    if (out) *out = nullptr;
    return protect([&] { if (!out || !path || !*path || overwrite > 1u) invalid_argument("BED sink arguments are invalid"); *out = new BedSinkHandle(path, n, overwrite != 0u, capacity); });
}
status_t bed_sink_append(handle_t* sink, const handle_t* h1, const handle_t* h2) {
    return protect([&] { required<BedSinkHandle>(sink, "BED sink").value.append(packed(h1), packed(h2)); });
}
status_t bed_sink_finalize(handle_t* sink) { return protect([&] { required<BedSinkHandle>(sink, "BED sink").value.finalize(); }); }
status_t bed_sink_info(const handle_t* sink, BedSinkInfo* out) {
    return protect([&] { if (!out) invalid_argument("BED sink info output is null"); const auto& x=required<BedSinkHandle>(sink,"BED sink").value; *out={x.individual_count(),x.variant_count(),x.bytes_written(),x.conversion_buffer_bytes(),x.lifecycle_object_bytes(),static_cast<int>(x.state())}; });
}
status_t bed_sink_close(handle_t* sink) { return protect([&] { delete &required<BedSinkHandle>(sink, "BED sink"); }); }

status_t hap_sink_create(const char* path, std::uint64_t n,
                         std::uint32_t overwrite, handle_t** out) {
    if (out) *out = nullptr;
    return protect([&] { if (!out || !path || !*path || overwrite > 1u) invalid_argument("HAP sink arguments are invalid"); *out = new HapSinkHandle(path, n, overwrite != 0u); });
}
status_t hap_sink_append(handle_t* sink, const handle_t* h1, const handle_t* h2) {
    return protect([&] { required<HapSinkHandle>(sink,"HAP sink").value.append(packed(h1),packed(h2)); });
}
status_t hap_sink_begin(handle_t* sink, std::uint64_t marker_count) {
    return protect([&] { required<HapSinkHandle>(sink,"HAP sink").value.begin_chromosome(marker_count); });
}
status_t hap_sink_write_batch(handle_t* sink, const handle_t* h1,
                              const handle_t* h2, std::uint64_t offset) {
    return protect([&] { required<HapSinkHandle>(sink,"HAP sink").value.write_batch(packed(h1),packed(h2),offset); });
}
status_t hap_sink_finalize(handle_t* sink) { return protect([&] { required<HapSinkHandle>(sink,"HAP sink").value.finalize(); }); }
status_t hap_sink_info(const handle_t* sink, HapSinkInfo* out) {
    return protect([&] { if(!out) invalid_argument("HAP sink info output is null"); const auto& x=required<HapSinkHandle>(sink,"HAP sink").value; *out={x.individual_count(),x.marker_count(),x.chromosome_count(),x.bytes_written(),static_cast<int>(x.state())}; });
}
status_t hap_sink_close(handle_t* sink) { return protect([&] { delete &required<HapSinkHandle>(sink,"HAP sink"); }); }
status_t hap_reader_open(const char* path, handle_t** out) {
    if(out) *out=nullptr;
    return protect([&] { if(!out || !path || !*path) invalid_argument("HAP reader path/output is empty"); *out=new HapReaderHandle(path); });
}
status_t hap_reader_close(handle_t* reader) { return protect([&] { delete &required<HapReaderHandle>(reader,"HAP reader"); }); }
status_t hap_reader_dimensions(const handle_t* reader, std::uint64_t* n,
                               std::uint64_t* m, std::uint64_t* c) {
    return protect([&] { if(!n||!m||!c) invalid_argument("HAP dimension output is null"); const auto& x=required<HapReaderHandle>(reader,"HAP reader").value; *n=x.individual_count(); *m=x.marker_count(); *c=x.chromosome_count(); });
}
status_t hap_reader_chromosome_info(const handle_t* reader, std::uint64_t index,
                                    HapChromosomeInfo* out) {
    return protect([&] { if(!out) invalid_argument("HAP chromosome output is null"); const auto& x=required<HapReaderHandle>(reader,"HAP reader").value.chromosome(index); *out={x.global_start_marker,x.marker_count,x.h1_offset,x.h2_offset,x.bytes_per_phase}; });
}
status_t hap_reader_load(const handle_t* reader, std::uint64_t index,
                         handle_t** h1, handle_t** h2) {
    if (h1) *h1 = nullptr;
    if (h2) *h2 = nullptr;
    return protect([&] { if(!h1||!h2||h1==h2) invalid_argument("HAP phase outputs are invalid"); auto x=required<HapReaderHandle>(reader,"HAP reader").value.load_chromosome(index); std::unique_ptr<PackedHandle> first(new PackedHandle(std::move(x.first))); std::unique_ptr<PackedHandle> second(new PackedHandle(std::move(x.second))); *h1=first.release(); *h2=second.release(); });
}

} // namespace gsim::native::api
