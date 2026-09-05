#ifndef GSIM_STANDALONE_BED_HPP
#define GSIM_STANDALONE_BED_HPP

#include "standalone_packed.h"

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

namespace gsim::native {

enum class BedSinkState : std::uint32_t {
    open = 0u,
    finalized = 1u,
    failed = 2u
};

class PhasedBedSink final {
public:
    PhasedBedSink(std::string destination_utf8,
                  std::uint64_t individual_count,
                  bool overwrite,
                  std::uint64_t buffer_variant_capacity);
    ~PhasedBedSink() noexcept;

    PhasedBedSink(const PhasedBedSink&) = delete;
    PhasedBedSink& operator=(const PhasedBedSink&) = delete;
    PhasedBedSink(PhasedBedSink&&) = delete;
    PhasedBedSink& operator=(PhasedBedSink&&) = delete;

    void append(const PhasedHaplotypeMatrix& h1,
                const PhasedHaplotypeMatrix& h2);
    void finalize();

    std::uint64_t individual_count() const noexcept {
        return individual_count_;
    }
    std::uint64_t variant_count() const noexcept { return variant_count_; }
    std::uint64_t bytes_written() const noexcept { return bytes_written_; }
    std::uint64_t conversion_buffer_bytes() const noexcept {
        return static_cast<std::uint64_t>(buffer_.size());
    }
    std::uint64_t lifecycle_object_bytes() const noexcept {
        return static_cast<std::uint64_t>(sizeof(*this));
    }
    BedSinkState state() const noexcept { return state_; }

private:
    void require_open(const char* operation) const;
    void create_temporary_file();
    void write_bytes(const std::uint8_t* values, std::size_t count,
                     const char* operation);
    void publish();
    void fail_and_cleanup() noexcept;
    std::string display(const std::filesystem::path& path) const;

    std::filesystem::path destination_;
    std::filesystem::path temporary_;
    std::FILE* stream_;
    std::uint64_t individual_count_;
    std::uint64_t bytes_per_variant_;
    std::uint64_t buffer_variant_capacity_;
    std::uint64_t variant_count_;
    std::uint64_t bytes_written_;
    bool overwrite_;
    BedSinkState state_;
    std::vector<std::uint8_t> buffer_;
};

} // namespace gsim::native

#endif
