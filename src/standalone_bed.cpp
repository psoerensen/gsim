#include "standalone_bed.h"

#include "standalone_error.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <limits>
#include <sstream>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <process.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace gsim::native {
namespace {

constexpr std::uint64_t header_size = 3u;
constexpr std::uint64_t maximum_buffer_variants = 4096u;
std::atomic<std::uint64_t> temporary_counter{0u};

std::uint64_t process_id() noexcept {
#ifdef _WIN32
    return static_cast<std::uint64_t>(_getpid());
#else
    return static_cast<std::uint64_t>(getpid());
#endif
}

int exclusive_open(const std::filesystem::path& path) noexcept {
#ifdef _WIN32
    return _wopen(path.c_str(), _O_BINARY | _O_WRONLY | _O_CREAT | _O_EXCL,
                  _S_IREAD | _S_IWRITE);
#else
    return ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL,
                  S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
#endif
}

void close_descriptor(int descriptor) noexcept {
#ifdef _WIN32
    (void)_close(descriptor);
#else
    (void)::close(descriptor);
#endif
}

std::FILE* descriptor_stream(int descriptor) noexcept {
#ifdef _WIN32
    return _fdopen(descriptor, "wb");
#else
    return fdopen(descriptor, "wb");
#endif
}

std::string errno_message() {
    return std::strerror(errno);
}

void remove_owned(const std::filesystem::path& path) noexcept {
    if (path.empty()) return;
    std::error_code ignored;
    (void)std::filesystem::remove(path, ignored);
}

} // namespace

PhasedBedSink::PhasedBedSink(std::string destination_utf8,
                             std::uint64_t individual_count,
                             bool overwrite,
                             std::uint64_t buffer_variant_capacity)
    : stream_(nullptr),
      individual_count_(individual_count),
      bytes_per_variant_(0u),
      buffer_variant_capacity_(buffer_variant_capacity),
      variant_count_(0u),
      bytes_written_(0u),
      overwrite_(overwrite),
      state_(BedSinkState::failed) {
    if (destination_utf8.empty()) {
        invalid_argument("BED sink destination path must not be empty");
    }
    if (individual_count_ == 0u) {
        invalid_argument("BED sink individual count must be positive");
    }
    if (buffer_variant_capacity_ == 0u ||
        buffer_variant_capacity_ > maximum_buffer_variants) {
        invalid_argument("BED sink buffer capacity must be in [1, 4096]");
    }

    std::error_code path_error;
    destination_ = std::filesystem::absolute(
        std::filesystem::u8path(destination_utf8), path_error);
    if (path_error || destination_.filename().empty()) {
        throw Error(GBITS_STATUS_INVALID_ARGUMENT,
                    "cannot resolve BED destination '" + destination_utf8 +
                        "': " + path_error.message());
    }
    destination_ = destination_.lexically_normal();
    const std::filesystem::path parent = destination_.parent_path();
    if (!std::filesystem::is_directory(parent, path_error) || path_error) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    "BED destination parent is not an accessible directory for '" +
                        display(destination_) + "'");
    }
    const bool destination_exists = std::filesystem::exists(destination_, path_error);
    if (path_error) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    "cannot inspect BED destination '" + display(destination_) +
                        "': " + path_error.message());
    }
    if (destination_exists && !overwrite_) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    "BED destination already exists and overwrite is disabled: '" +
                        display(destination_) + "'");
    }
    if (destination_exists) {
        const bool destination_is_directory =
            std::filesystem::is_directory(destination_, path_error);
        if (path_error) {
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot inspect BED destination type '" +
                            display(destination_) + "': " +
                            path_error.message());
        }
        if (destination_is_directory) {
            throw Error(GBITS_STATUS_IO_ERROR,
                        "BED destination is a directory: '" +
                            display(destination_) + "'");
        }
    }

    bytes_per_variant_ = individual_count_ / 4u +
                         static_cast<std::uint64_t>(individual_count_ % 4u != 0u);
    const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
    if (buffer_variant_capacity_ > maximum / bytes_per_variant_) {
        invalid_argument("BED conversion-buffer size overflows uint64_t");
    }
    const std::uint64_t buffer_bytes =
        buffer_variant_capacity_ * bytes_per_variant_;
    if (buffer_bytes > static_cast<std::uint64_t>(
                           std::numeric_limits<std::size_t>::max())) {
        invalid_argument("BED conversion buffer exceeds platform allocation limits");
    }
    buffer_.assign(static_cast<std::size_t>(buffer_bytes), 0u);

    try {
        create_temporary_file();
        const std::uint8_t header[3] = {0x6cu, 0x1bu, 0x01u};
        write_bytes(header, sizeof(header), "write BED header");
        bytes_written_ = header_size;
        state_ = BedSinkState::open;
    } catch (...) {
        fail_and_cleanup();
        throw;
    }
}

PhasedBedSink::~PhasedBedSink() noexcept {
    if (state_ != BedSinkState::finalized) fail_and_cleanup();
}

std::string PhasedBedSink::display(const std::filesystem::path& path) const {
    return path.u8string();
}

void PhasedBedSink::create_temporary_file() {
    for (unsigned int attempt = 0u; attempt < 128u; ++attempt) {
        const std::uint64_t serial =
            temporary_counter.fetch_add(1u, std::memory_order_relaxed);
        std::ostringstream suffix;
        suffix << destination_.filename().u8string() << ".gsim.tmp."
               << process_id() << '.' << serial;
        temporary_ = destination_.parent_path() /
                     std::filesystem::u8path(suffix.str());
        errno = 0;
        const int descriptor = exclusive_open(temporary_);
        if (descriptor < 0) {
            if (errno == EEXIST) continue;
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot create BED temporary file '" +
                            display(temporary_) + "': " + errno_message());
        }
        stream_ = descriptor_stream(descriptor);
        if (stream_ == nullptr) {
            const std::string detail = errno_message();
            close_descriptor(descriptor);
            remove_owned(temporary_);
            temporary_.clear();
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot open BED temporary stream for '" +
                            display(destination_) + "': " + detail);
        }
        return;
    }
    throw Error(GBITS_STATUS_IO_ERROR,
                "cannot acquire an exclusive BED temporary file for '" +
                    display(destination_) + "'");
}

void PhasedBedSink::write_bytes(const std::uint8_t* values, std::size_t count,
                                const char* operation) {
    if (count == 0u) return;
    if (stream_ == nullptr ||
        std::fwrite(values, 1u, count, stream_) != count) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    std::string(operation) + " failed for BED destination '" +
                        display(destination_) + "'");
    }
}

void PhasedBedSink::require_open(const char* operation) const {
    if (state_ == BedSinkState::finalized) {
        throw Error(GBITS_STATUS_INVALID_ARGUMENT,
                    std::string(operation) +
                        " is invalid after BED sink finalization");
    }
    if (state_ == BedSinkState::failed) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    std::string(operation) + " is invalid for a failed BED sink");
    }
}

void PhasedBedSink::append(const PhasedHaplotypeMatrix& h1,
                           const PhasedHaplotypeMatrix& h2) {
    require_open("BED append");
    try {
        if (h1.individual_count() != h2.individual_count() ||
            h1.marker_count() != h2.marker_count() ||
            h1.words_per_marker() != h2.words_per_marker()) {
            invalid_argument("BED append requires compatible H1/H2 dimensions and stride");
        }
        if (h1.individual_count() != individual_count_) {
            invalid_argument("BED append individual count differs from sink contract");
        }
        if (!h1.has_canonical_padding() || !h2.has_canonical_padding()) {
            invalid_argument("BED append requires canonical zero haplotype padding");
        }

        const std::uint64_t markers = h1.marker_count();
        const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
        if (markers > maximum - variant_count_ ||
            markers > (maximum - bytes_written_) / bytes_per_variant_) {
            invalid_argument("BED append dimensions overflow counters");
        }
        const std::uint64_t resulting_bytes =
            bytes_written_ + markers * bytes_per_variant_;
        if (resulting_bytes > static_cast<std::uint64_t>(
                                  std::numeric_limits<std::streamoff>::max())) {
            invalid_argument("BED output exceeds platform file-offset limits");
        }

        for (std::uint64_t first = 0u; first < markers;) {
            const std::uint64_t count =
                std::min(buffer_variant_capacity_, markers - first);
            const std::uint64_t block_bytes = count * bytes_per_variant_;
            std::fill_n(buffer_.data(),
                        static_cast<std::size_t>(block_bytes),
                        std::uint8_t{0});
            for (std::uint64_t local = 0u; local < count; ++local) {
                const std::uint64_t marker = first + local;
                std::uint8_t* record = buffer_.data() +
                    static_cast<std::size_t>(local * bytes_per_variant_);
                for (std::uint64_t individual = 0u;
                     individual < individual_count_; ++individual) {
                    const std::uint8_t dosage = static_cast<std::uint8_t>(
                        h1.allele(individual, marker) +
                        h2.allele(individual, marker));
                    const std::uint8_t code =
                        dosage == 0u ? std::uint8_t{0x03u}
                        : dosage == 1u ? std::uint8_t{0x02u}
                                        : std::uint8_t{0x00u};
                    const unsigned int shift =
                        static_cast<unsigned int>(2u * (individual % 4u));
                    record[static_cast<std::size_t>(individual / 4u)] |=
                        static_cast<std::uint8_t>(code << shift);
                }
            }
            write_bytes(buffer_.data(), static_cast<std::size_t>(block_bytes),
                        "append packed chromosome");
            first += count;
        }
        variant_count_ += markers;
        bytes_written_ = resulting_bytes;
    } catch (...) {
        fail_and_cleanup();
        throw;
    }
}

void PhasedBedSink::publish() {
#ifdef _WIN32
    DWORD flags = MOVEFILE_WRITE_THROUGH;
    if (overwrite_) flags |= MOVEFILE_REPLACE_EXISTING;
    if (!MoveFileExW(temporary_.c_str(), destination_.c_str(), flags)) {
        std::ostringstream message;
        message << "cannot publish BED destination '" << display(destination_)
                << "' (Windows error " << GetLastError() << ')';
        throw Error(GBITS_STATUS_IO_ERROR, message.str());
    }
#else
    if (overwrite_) {
        if (::rename(temporary_.c_str(), destination_.c_str()) != 0) {
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot replace BED destination '" +
                            display(destination_) + "': " + errno_message());
        }
    } else {
        if (::link(temporary_.c_str(), destination_.c_str()) != 0) {
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot publish BED destination '" +
                            display(destination_) + "': " + errno_message());
        }
        if (::unlink(temporary_.c_str()) != 0) {
            std::error_code rollback_error;
            (void)std::filesystem::remove(destination_, rollback_error);
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot remove BED temporary link after publishing '" +
                            display(destination_) + "': " + errno_message());
        }
    }
#endif
}

void PhasedBedSink::finalize() {
    require_open("BED finalize");
    try {
        if (stream_ == nullptr || std::fflush(stream_) != 0) {
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot flush BED temporary file for '" +
                            display(destination_) + "'");
        }
        const int close_status = std::fclose(stream_);
        stream_ = nullptr;
        if (close_status != 0) {
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot close BED temporary file for '" +
                            display(destination_) + "'");
        }
        publish();
        temporary_.clear();
        state_ = BedSinkState::finalized;
    } catch (...) {
        fail_and_cleanup();
        throw;
    }
}

void PhasedBedSink::fail_and_cleanup() noexcept {
    if (stream_ != nullptr) {
        (void)std::fclose(stream_);
        stream_ = nullptr;
    }
    remove_owned(temporary_);
    temporary_.clear();
    if (state_ != BedSinkState::finalized) state_ = BedSinkState::failed;
}

} // namespace gsim::native
