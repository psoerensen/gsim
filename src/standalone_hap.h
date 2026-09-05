#ifndef GSIM_STANDALONE_HAP_HPP
#define GSIM_STANDALONE_HAP_HPP

#include "standalone_packed.h"

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace gsim::native {

enum class HapSinkState : std::uint32_t {
    open = 0u,
    finalized = 1u,
    failed = 2u
};

struct HapChromosomeInfo final {
    std::uint64_t global_start_marker;
    std::uint64_t marker_count;
    std::uint64_t h1_offset;
    std::uint64_t h2_offset;
    std::uint64_t bytes_per_phase;
};

class PhasedHapSink final {
public:
    PhasedHapSink(std::string destination_utf8,
                  std::uint64_t individual_count, bool overwrite);
    ~PhasedHapSink() noexcept;

    PhasedHapSink(const PhasedHapSink&) = delete;
    PhasedHapSink& operator=(const PhasedHapSink&) = delete;

    void append(const PhasedHaplotypeMatrix& h1,
                const PhasedHaplotypeMatrix& h2);
    void finalize();

    std::uint64_t individual_count() const noexcept { return individual_count_; }
    std::uint64_t marker_count() const noexcept { return marker_count_; }
    std::uint64_t chromosome_count() const noexcept {
        return static_cast<std::uint64_t>(chromosomes_.size());
    }
    std::uint64_t bytes_written() const noexcept { return bytes_written_; }
    HapSinkState state() const noexcept { return state_; }

private:
    void require_open(const char* operation) const;
    void create_temporary_file();
    void write_bytes(const std::uint8_t* values, std::size_t count,
                     const char* operation);
    void write_u32(std::uint32_t value, const char* operation);
    void write_u64(std::uint64_t value, const char* operation);
    void seek(std::uint64_t offset, const char* operation);
    void publish();
    void fail_and_cleanup() noexcept;
    std::string display(const std::filesystem::path& path) const;

    std::filesystem::path destination_;
    std::filesystem::path temporary_;
    std::FILE* stream_;
    std::uint64_t individual_count_;
    std::uint64_t marker_count_;
    std::uint64_t bytes_written_;
    bool overwrite_;
    HapSinkState state_;
    std::vector<HapChromosomeInfo> chromosomes_;
};

class PhasedHapReader final {
public:
    explicit PhasedHapReader(std::string path_utf8);

    std::uint64_t individual_count() const noexcept { return individual_count_; }
    std::uint64_t marker_count() const noexcept { return marker_count_; }
    std::uint64_t chromosome_count() const noexcept {
        return static_cast<std::uint64_t>(chromosomes_.size());
    }
    const HapChromosomeInfo& chromosome(std::uint64_t index) const;
    std::pair<PhasedHaplotypeMatrix, PhasedHaplotypeMatrix>
    load_chromosome(std::uint64_t index) const;

private:
    std::filesystem::path path_;
    std::uint64_t individual_count_;
    std::uint64_t marker_count_;
    std::uint64_t file_size_;
    std::vector<HapChromosomeInfo> chromosomes_;
};

} // namespace gsim::native

#endif
