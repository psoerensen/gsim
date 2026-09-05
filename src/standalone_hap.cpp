#include "standalone_hap.h"

#include "standalone_error.h"

#include <array>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <fstream>
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

constexpr std::uint64_t header_size = 64u;
constexpr std::uint64_t table_entry_size = 48u;
constexpr std::uint64_t data_offset = header_size;
std::atomic<std::uint64_t> temporary_counter{0u};

std::uint64_t checked_add(std::uint64_t left, std::uint64_t right,
                          const char* name) {
    if (right > std::numeric_limits<std::uint64_t>::max() - left) {
        invalid_argument(std::string(name) + " overflows uint64_t");
    }
    return left + right;
}

std::uint64_t checked_multiply(std::uint64_t left, std::uint64_t right,
                               const char* name) {
    if (left != 0u && right > std::numeric_limits<std::uint64_t>::max() / left) {
        invalid_argument(std::string(name) + " overflows uint64_t");
    }
    return left * right;
}

std::array<std::uint8_t, 4> encode_u32(std::uint32_t value) {
    return {{static_cast<std::uint8_t>(value),
             static_cast<std::uint8_t>(value >> 8u),
             static_cast<std::uint8_t>(value >> 16u),
             static_cast<std::uint8_t>(value >> 24u)}};
}

std::array<std::uint8_t, 8> encode_u64(std::uint64_t value) {
    std::array<std::uint8_t, 8> bytes{};
    for (unsigned int i = 0u; i < 8u; ++i) {
        bytes[i] = static_cast<std::uint8_t>(value >> (8u * i));
    }
    return bytes;
}

std::uint32_t decode_u32(const std::uint8_t* bytes) {
    std::uint32_t value = 0u;
    for (unsigned int i = 0u; i < 4u; ++i) {
        value |= static_cast<std::uint32_t>(bytes[i]) << (8u * i);
    }
    return value;
}

std::uint64_t decode_u64(const std::uint8_t* bytes) {
    std::uint64_t value = 0u;
    for (unsigned int i = 0u; i < 8u; ++i) {
        value |= static_cast<std::uint64_t>(bytes[i]) << (8u * i);
    }
    return value;
}

std::uint64_t process_id() noexcept {
#ifdef _WIN32
    return static_cast<std::uint64_t>(_getpid());
#else
    return static_cast<std::uint64_t>(getpid());
#endif
}

int exclusive_open(const std::filesystem::path& path) noexcept {
#ifdef _WIN32
    return _wopen(path.c_str(), _O_BINARY | _O_RDWR | _O_CREAT | _O_EXCL,
                  _S_IREAD | _S_IWRITE);
#else
    return ::open(path.c_str(), O_RDWR | O_CREAT | O_EXCL,
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
    return _fdopen(descriptor, "w+b");
#else
    return fdopen(descriptor, "w+b");
#endif
}

int seek_file(std::FILE* stream, std::uint64_t offset) noexcept {
#ifdef _WIN32
    if (offset > static_cast<std::uint64_t>(
                     std::numeric_limits<__int64>::max())) return -1;
    return _fseeki64(stream, static_cast<__int64>(offset), SEEK_SET);
#else
    if (offset > static_cast<std::uint64_t>(
                     std::numeric_limits<off_t>::max())) return -1;
    return fseeko(stream, static_cast<off_t>(offset), SEEK_SET);
#endif
}

std::string errno_message() { return std::strerror(errno); }

void remove_owned(const std::filesystem::path& path) noexcept {
    if (path.empty()) return;
    std::error_code ignored;
    (void)std::filesystem::remove(path, ignored);
}

void read_exact(std::ifstream& input, std::uint64_t offset,
                std::uint8_t* output, std::size_t count,
                const std::string& context) {
    if (offset > static_cast<std::uint64_t>(
                     std::numeric_limits<std::streamoff>::max())) {
        invalid_argument("HAP file offset exceeds platform stream limits");
    }
    input.clear();
    input.seekg(static_cast<std::streamoff>(offset), std::ios::beg);
    if (!input || !input.read(reinterpret_cast<char*>(output),
                              static_cast<std::streamsize>(count))) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    "cannot read " + context + " from HAP file");
    }
}

} // namespace

PhasedHapSink::PhasedHapSink(std::string destination_utf8,
                             std::uint64_t individual_count, bool overwrite)
    : stream_(nullptr), individual_count_(individual_count), marker_count_(0u),
      bytes_written_(0u), overwrite_(overwrite), state_(HapSinkState::failed) {
    if (destination_utf8.empty()) invalid_argument("HAP destination path must not be empty");
    if (individual_count_ == 0u) invalid_argument("HAP sample count must be positive");
    std::error_code error;
    destination_ = std::filesystem::absolute(
        std::filesystem::u8path(destination_utf8), error);
    if (error || destination_.filename().empty()) {
        throw Error(GBITS_STATUS_INVALID_ARGUMENT,
                    "cannot resolve HAP destination '" + destination_utf8 + "'");
    }
    destination_ = destination_.lexically_normal();
    if (!std::filesystem::is_directory(destination_.parent_path(), error) || error) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    "HAP destination parent is not an accessible directory for '" +
                        display(destination_) + "'");
    }
    const bool exists = std::filesystem::exists(destination_, error);
    if (error) throw Error(GBITS_STATUS_IO_ERROR, "cannot inspect HAP destination");
    if (exists && !overwrite_) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    "HAP destination already exists and overwrite is disabled: '" +
                        display(destination_) + "'");
    }
    if (exists && std::filesystem::is_directory(destination_, error)) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    "HAP destination is a directory: '" + display(destination_) + "'");
    }
    try {
        create_temporary_file();
        std::array<std::uint8_t, static_cast<std::size_t>(header_size)> empty{};
        write_bytes(empty.data(), empty.size(), "write HAP placeholder header");
        bytes_written_ = header_size;
        state_ = HapSinkState::open;
    } catch (...) {
        fail_and_cleanup();
        throw;
    }
}

PhasedHapSink::~PhasedHapSink() noexcept {
    if (state_ != HapSinkState::finalized) fail_and_cleanup();
}

std::string PhasedHapSink::display(const std::filesystem::path& path) const {
    return path.u8string();
}

void PhasedHapSink::create_temporary_file() {
    for (unsigned int attempt = 0u; attempt < 128u; ++attempt) {
        const std::uint64_t serial = temporary_counter.fetch_add(1u);
        std::ostringstream name;
        name << destination_.filename().u8string() << ".gsim.hap.tmp."
             << process_id() << '.' << serial;
        temporary_ = destination_.parent_path() / std::filesystem::u8path(name.str());
        errno = 0;
        const int descriptor = exclusive_open(temporary_);
        if (descriptor < 0) {
            if (errno == EEXIST) continue;
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot create HAP temporary file '" + display(temporary_) +
                            "': " + errno_message());
        }
        stream_ = descriptor_stream(descriptor);
        if (stream_ == nullptr) {
            const std::string detail = errno_message();
            close_descriptor(descriptor);
            remove_owned(temporary_);
            temporary_.clear();
            throw Error(GBITS_STATUS_IO_ERROR,
                        "cannot open HAP temporary stream: " + detail);
        }
        return;
    }
    throw Error(GBITS_STATUS_IO_ERROR, "cannot acquire an exclusive HAP temporary file");
}

void PhasedHapSink::write_bytes(const std::uint8_t* values, std::size_t count,
                                const char* operation) {
    if (count != 0u && (stream_ == nullptr ||
        std::fwrite(values, 1u, count, stream_) != count)) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    std::string(operation) + " failed for HAP destination '" +
                        display(destination_) + "'");
    }
}

void PhasedHapSink::write_u32(std::uint32_t value, const char* operation) {
    const auto bytes = encode_u32(value);
    write_bytes(bytes.data(), bytes.size(), operation);
}

void PhasedHapSink::write_u64(std::uint64_t value, const char* operation) {
    const auto bytes = encode_u64(value);
    write_bytes(bytes.data(), bytes.size(), operation);
}

void PhasedHapSink::seek(std::uint64_t offset, const char* operation) {
    if (stream_ == nullptr || seek_file(stream_, offset) != 0) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    std::string(operation) + " failed for HAP destination '" +
                        display(destination_) + "'");
    }
}

void PhasedHapSink::require_open(const char* operation) const {
    if (state_ == HapSinkState::finalized) {
        invalid_argument(std::string(operation) + " is invalid after HAP finalization");
    }
    if (state_ == HapSinkState::failed) {
        throw Error(GBITS_STATUS_IO_ERROR,
                    std::string(operation) + " is invalid for a failed HAP sink");
    }
}

void PhasedHapSink::append(const PhasedHaplotypeMatrix& h1,
                           const PhasedHaplotypeMatrix& h2) {
    require_open("HAP append");
    try {
        if (h1.individual_count() != h2.individual_count() ||
            h1.marker_count() != h2.marker_count() ||
            h1.words_per_marker() != h2.words_per_marker()) {
            invalid_argument("HAP append requires compatible H1/H2 dimensions and stride");
        }
        if (h1.individual_count() != individual_count_) {
            invalid_argument("HAP append sample count differs from sink contract");
        }
        if (!h1.has_canonical_padding() || !h2.has_canonical_padding()) {
            invalid_argument("HAP append requires canonical zero padding");
        }
        const std::uint64_t phase_bytes = h1.storage_bytes();
        const std::uint64_t next_markers = checked_add(
            marker_count_, h1.marker_count(), "HAP marker count");
        const std::uint64_t h1_offset = bytes_written_;
        const std::uint64_t h2_offset = checked_add(h1_offset, phase_bytes,
                                                    "HAP H2 offset");
        const std::uint64_t next_bytes = checked_add(
            h2_offset, phase_bytes, "HAP data size");
        if (next_bytes > static_cast<std::uint64_t>(
                             std::numeric_limits<std::streamoff>::max())) {
            invalid_argument("HAP data exceeds platform file-offset limits");
        }
        HapChromosomeInfo info{marker_count_, h1.marker_count(), h1_offset,
                               h2_offset, phase_bytes};
        for (const PhasedHaplotypeMatrix* phase : {&h1, &h2}) {
            for (std::uint64_t marker = 0u; marker < phase->marker_count(); ++marker) {
                for (std::uint64_t word = 0u; word < phase->words_per_marker(); ++word) {
                    write_u64(phase->word(marker, word), "append packed HAP word");
                }
            }
        }
        chromosomes_.push_back(info);
        marker_count_ = next_markers;
        bytes_written_ = next_bytes;
    } catch (...) {
        fail_and_cleanup();
        throw;
    }
}

void PhasedHapSink::publish() {
#ifdef _WIN32
    DWORD flags = MOVEFILE_WRITE_THROUGH;
    if (overwrite_) flags |= MOVEFILE_REPLACE_EXISTING;
    if (!MoveFileExW(temporary_.c_str(), destination_.c_str(), flags)) {
        std::ostringstream message;
        message << "cannot publish HAP destination '" << display(destination_)
                << "' (Windows error " << GetLastError() << ')';
        throw Error(GBITS_STATUS_IO_ERROR, message.str());
    }
#else
    if (overwrite_) {
        if (::rename(temporary_.c_str(), destination_.c_str()) != 0) {
            throw Error(GBITS_STATUS_IO_ERROR, "cannot replace HAP destination: " + errno_message());
        }
    } else {
        if (::link(temporary_.c_str(), destination_.c_str()) != 0) {
            throw Error(GBITS_STATUS_IO_ERROR, "cannot publish HAP destination: " + errno_message());
        }
        if (::unlink(temporary_.c_str()) != 0) {
            std::error_code ignored;
            (void)std::filesystem::remove(destination_, ignored);
            throw Error(GBITS_STATUS_IO_ERROR, "cannot unlink HAP temporary after publish");
        }
    }
#endif
}

void PhasedHapSink::finalize() {
    require_open("HAP finalize");
    try {
        if (chromosomes_.empty() || marker_count_ == 0u) {
            invalid_argument("HAP finalization requires at least one nonempty chromosome");
        }
        const std::uint64_t table_offset = bytes_written_;
        const std::uint64_t table_bytes = checked_multiply(
            static_cast<std::uint64_t>(chromosomes_.size()), table_entry_size,
            "HAP chromosome table size");
        const std::uint64_t final_size = checked_add(table_offset, table_bytes,
                                                     "HAP file size");
        if (final_size > static_cast<std::uint64_t>(
                             std::numeric_limits<std::streamoff>::max())) {
            invalid_argument("HAP file exceeds platform file-offset limits");
        }
        for (const HapChromosomeInfo& entry : chromosomes_) {
            write_u64(entry.global_start_marker, "write HAP table");
            write_u64(entry.marker_count, "write HAP table");
            write_u64(entry.h1_offset, "write HAP table");
            write_u64(entry.h2_offset, "write HAP table");
            write_u64(entry.bytes_per_phase, "write HAP table");
            write_u64(0u, "write HAP table reserved field");
        }
        seek(0u, "seek HAP header");
        const std::uint8_t signature[4] = {0x48u, 0x41u, 0x50u, 0x01u};
        write_bytes(signature, 4u, "write HAP signature");
        write_u32(static_cast<std::uint32_t>(header_size), "write HAP header size");
        write_u32(0u, "write HAP flags");
        write_u32(0u, "write HAP reserved field");
        write_u64(individual_count_, "write HAP sample count");
        write_u64(marker_count_, "write HAP marker count");
        write_u64(static_cast<std::uint64_t>(chromosomes_.size()),
                  "write HAP chromosome count");
        write_u64(table_offset, "write HAP table offset");
        write_u64(data_offset, "write HAP data offset");
        write_u64(0u, "write HAP reserved field");
        if (std::fflush(stream_) != 0) {
            throw Error(GBITS_STATUS_IO_ERROR, "cannot flush HAP temporary file");
        }
        if (std::fclose(stream_) != 0) {
            stream_ = nullptr;
            throw Error(GBITS_STATUS_IO_ERROR, "cannot close HAP temporary file");
        }
        stream_ = nullptr;
        std::error_code size_error;
        const std::uint64_t observed = std::filesystem::file_size(temporary_, size_error);
        if (size_error || observed != final_size) {
            throw Error(GBITS_STATUS_IO_ERROR, "HAP temporary file size verification failed");
        }
        bytes_written_ = final_size;
        publish();
        temporary_.clear();
        state_ = HapSinkState::finalized;
    } catch (...) {
        fail_and_cleanup();
        throw;
    }
}

void PhasedHapSink::fail_and_cleanup() noexcept {
    if (stream_ != nullptr) {
        (void)std::fclose(stream_);
        stream_ = nullptr;
    }
    remove_owned(temporary_);
    temporary_.clear();
    if (state_ != HapSinkState::finalized) state_ = HapSinkState::failed;
}

PhasedHapReader::PhasedHapReader(std::string path_utf8)
    : individual_count_(0u), marker_count_(0u), file_size_(0u) {
    if (path_utf8.empty()) invalid_argument("HAP input path must not be empty");
    std::error_code error;
    path_ = std::filesystem::absolute(std::filesystem::u8path(path_utf8), error);
    if (error) throw Error(GBITS_STATUS_INVALID_ARGUMENT, "cannot resolve HAP input path");
    file_size_ = std::filesystem::file_size(path_, error);
    if (error) throw Error(GBITS_STATUS_IO_ERROR, "cannot inspect HAP input file");
    std::ifstream input(path_, std::ios::binary);
    if (!input) throw Error(GBITS_STATUS_IO_ERROR, "cannot open HAP input file");
    std::array<std::uint8_t, static_cast<std::size_t>(header_size)> header{};
    if (file_size_ < 4u) throw Error(GBITS_STATUS_IO_ERROR, "HAP file is truncated before signature");
    read_exact(input, 0u, header.data(), 4u, "HAP signature");
    if (header[0] != 0x48u || header[1] != 0x41u || header[2] != 0x50u) {
        invalid_argument("HAP file signature is not HAP");
    }
    if (header[3] != 0x01u) invalid_argument("unsupported HAP format version");
    if (file_size_ < header_size) throw Error(GBITS_STATUS_IO_ERROR, "HAP header is truncated");
    read_exact(input, 4u, header.data() + 4u,
               static_cast<std::size_t>(header_size - 4u), "HAP header");
    if (decode_u32(header.data() + 4u) != header_size ||
        decode_u32(header.data() + 8u) != 0u ||
        decode_u32(header.data() + 12u) != 0u ||
        decode_u64(header.data() + 56u) != 0u) {
        invalid_argument("HAP header has incompatible size, flags, or reserved fields");
    }
    individual_count_ = decode_u64(header.data() + 16u);
    marker_count_ = decode_u64(header.data() + 24u);
    const std::uint64_t chromosome_count = decode_u64(header.data() + 32u);
    const std::uint64_t table_offset = decode_u64(header.data() + 40u);
    if (decode_u64(header.data() + 48u) != data_offset || individual_count_ == 0u ||
        marker_count_ == 0u || chromosome_count == 0u) {
        invalid_argument("HAP header has invalid dimensions or data offset");
    }
    const std::uint64_t table_bytes = checked_multiply(
        chromosome_count, table_entry_size, "HAP chromosome table size");
    if (checked_add(table_offset, table_bytes, "HAP file size") != file_size_) {
        invalid_argument("HAP file length does not match its chromosome table");
    }
    if (chromosome_count > static_cast<std::uint64_t>(
                               std::numeric_limits<std::size_t>::max() /
                               sizeof(HapChromosomeInfo))) {
        invalid_argument("HAP chromosome table exceeds platform limits");
    }
    chromosomes_.reserve(static_cast<std::size_t>(chromosome_count));
    std::uint64_t expected_marker = 0u;
    std::uint64_t expected_offset = data_offset;
    std::array<std::uint8_t, static_cast<std::size_t>(table_entry_size)> entry{};
    const std::uint64_t words_per_marker = individual_count_ / 64u +
        static_cast<std::uint64_t>(individual_count_ % 64u != 0u);
    const std::uint64_t marker_stride = checked_multiply(
        words_per_marker, 8u, "HAP marker stride");
    for (std::uint64_t index = 0u; index < chromosome_count; ++index) {
        read_exact(input, checked_add(table_offset,
                   checked_multiply(index, table_entry_size, "HAP table index"),
                   "HAP table offset"), entry.data(), entry.size(), "HAP chromosome table");
        HapChromosomeInfo info{decode_u64(entry.data()), decode_u64(entry.data() + 8u),
                               decode_u64(entry.data() + 16u), decode_u64(entry.data() + 24u),
                               decode_u64(entry.data() + 32u)};
        if (decode_u64(entry.data() + 40u) != 0u || info.marker_count == 0u ||
            info.global_start_marker != expected_marker || info.h1_offset != expected_offset) {
            invalid_argument("HAP chromosome table has invalid range, offset, or reserved field");
        }
        const std::uint64_t expected_phase = checked_multiply(
            checked_multiply(info.marker_count, words_per_marker,
                             "HAP chromosome word count"), 8u,
            "HAP chromosome byte count");
        if (info.bytes_per_phase != expected_phase ||
            info.h2_offset != checked_add(info.h1_offset, expected_phase, "HAP H2 offset")) {
            invalid_argument("HAP chromosome table has inconsistent phase layout");
        }
        expected_marker = checked_add(expected_marker, info.marker_count,
                                      "HAP marker ranges");
        expected_offset = checked_add(info.h2_offset, expected_phase,
                                      "HAP chromosome data end");
        chromosomes_.push_back(info);
    }
    if (expected_marker != marker_count_ || expected_offset != table_offset) {
        invalid_argument("HAP chromosome ranges do not cover the declared data");
    }
    if (individual_count_ % 64u != 0u) {
        const unsigned int used = static_cast<unsigned int>(individual_count_ % 64u);
        const std::uint64_t padding_mask = ~((std::uint64_t{1} << used) - 1u);
        std::array<std::uint8_t, 8> word_bytes{};
        for (const HapChromosomeInfo& info : chromosomes_) {
            for (const std::uint64_t plane : {info.h1_offset, info.h2_offset}) {
                for (std::uint64_t marker = 0u; marker < info.marker_count; ++marker) {
                    const std::uint64_t last_word_offset = checked_add(
                        plane, checked_multiply(marker, marker_stride,
                                                "HAP padding offset"),
                        "HAP padding offset");
                    read_exact(input, checked_add(last_word_offset,
                        (words_per_marker - 1u) * 8u, "HAP padding word offset"),
                        word_bytes.data(), word_bytes.size(), "HAP padding word");
                    if ((decode_u64(word_bytes.data()) & padding_mask) != 0u) {
                        invalid_argument("HAP file contains noncanonical padding bits");
                    }
                }
            }
        }
    }
}

const HapChromosomeInfo& PhasedHapReader::chromosome(std::uint64_t index) const {
    if (index >= chromosomes_.size()) {
        throw Error(GBITS_STATUS_OUT_OF_RANGE, "HAP chromosome index is out of range");
    }
    return chromosomes_[static_cast<std::size_t>(index)];
}

std::pair<PhasedHaplotypeMatrix, PhasedHaplotypeMatrix>
PhasedHapReader::load_chromosome(std::uint64_t index) const {
    const HapChromosomeInfo& info = chromosome(index);
    PhasedHaplotypeMatrix h1(individual_count_, info.marker_count);
    PhasedHaplotypeMatrix h2(individual_count_, info.marker_count);
    std::ifstream input(path_, std::ios::binary);
    if (!input) throw Error(GBITS_STATUS_IO_ERROR, "cannot reopen HAP input file");
    std::array<std::uint8_t, 8> bytes{};
    for (auto phase : {std::make_pair(info.h1_offset, &h1),
                       std::make_pair(info.h2_offset, &h2)}) {
        std::uint64_t offset = phase.first;
        for (std::uint64_t marker = 0u; marker < info.marker_count; ++marker) {
            for (std::uint64_t word = 0u; word < phase.second->words_per_marker(); ++word) {
                read_exact(input, offset, bytes.data(), bytes.size(), "HAP packed word");
                phase.second->set_word(marker, word, decode_u64(bytes.data()));
                offset = checked_add(offset, 8u, "HAP packed word offset");
            }
        }
    }
    return std::make_pair(std::move(h1), std::move(h2));
}

} // namespace gsim::native
