#ifndef GSIM_STANDALONE_GBITS_API_HPP
#define GSIM_STANDALONE_GBITS_API_HPP

#include <cstdint>

namespace gsim::native::api {

using status_t = int;
using handle_t = void;

struct BedSinkInfo {
    std::uint64_t individual_count, variant_count, bytes_written;
    std::uint64_t conversion_buffer_bytes, lifecycle_object_bytes;
    int state;
};
struct HapSinkInfo {
    std::uint64_t individual_count, marker_count, chromosome_count, bytes_written;
    int state;
};
struct HapChromosomeInfo {
    std::uint64_t global_start_marker, marker_count, h1_offset, h2_offset;
    std::uint64_t bytes_per_phase;
};

std::uint32_t abi_version() noexcept;
const char* library_version() noexcept;
const char* last_error() noexcept;
status_t create_zero(std::uint64_t, std::uint64_t, handle_t**);
status_t create_values(std::uint64_t, std::uint64_t, const std::uint8_t*,
                       std::uint64_t, std::uint64_t, handle_t**);
status_t close(handle_t*);
status_t individual_count(const handle_t*, std::uint64_t*);
status_t marker_count(const handle_t*, std::uint64_t*);
status_t words_per_marker(const handle_t*, std::uint64_t*);
status_t storage_bytes(const handle_t*, std::uint64_t*);
status_t word(const handle_t*, std::uint64_t, std::uint64_t, std::uint64_t*);
status_t allele(const handle_t*, std::uint64_t, std::uint64_t, std::uint8_t*);
status_t set_allele(handle_t*, std::uint64_t, std::uint64_t, std::uint8_t);
status_t unpack(const handle_t*, std::uint8_t*, std::uint64_t, std::uint64_t);
status_t copy_interval(handle_t*, std::uint64_t, const handle_t*, std::uint64_t,
                       std::uint64_t, std::uint64_t);
status_t copy_filtered(handle_t*, std::uint64_t, const handle_t*, std::uint64_t,
                       std::uint64_t, std::uint64_t, double, const double*,
                       std::uint64_t);
status_t make_gamete(handle_t*, std::uint64_t, const handle_t*, const handle_t*,
                     std::uint64_t, std::uint32_t, const std::uint64_t*,
                     std::uint64_t);
status_t decode_genotypes(const handle_t*, const handle_t*, std::uint8_t*,
                          std::uint64_t, std::uint64_t);
status_t bed_open(const char*, std::uint64_t, std::uint64_t, handle_t**);
status_t bed_close(handle_t*);
status_t bed_read_variant(handle_t*, std::uint64_t, std::int8_t*, std::uint64_t);
status_t bed_sink_create(const char*, std::uint64_t, std::uint32_t,
                         std::uint64_t, handle_t**);
status_t bed_sink_append(handle_t*, const handle_t*, const handle_t*);
status_t bed_sink_finalize(handle_t*);
status_t bed_sink_info(const handle_t*, BedSinkInfo*);
status_t bed_sink_close(handle_t*);
status_t hap_sink_create(const char*, std::uint64_t, std::uint32_t, handle_t**);
status_t hap_sink_append(handle_t*, const handle_t*, const handle_t*);
status_t hap_sink_finalize(handle_t*);
status_t hap_sink_info(const handle_t*, HapSinkInfo*);
status_t hap_sink_close(handle_t*);
status_t hap_reader_open(const char*, handle_t**);
status_t hap_reader_close(handle_t*);
status_t hap_reader_dimensions(const handle_t*, std::uint64_t*, std::uint64_t*,
                               std::uint64_t*);
status_t hap_reader_chromosome_info(const handle_t*, std::uint64_t,
                                    HapChromosomeInfo*);
status_t hap_reader_load(const handle_t*, std::uint64_t, handle_t**, handle_t**);

} // namespace gsim::native::api

#endif
