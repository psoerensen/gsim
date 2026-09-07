// Adapted from gbits/src/bed_reader.cpp: checked record I/O and xmat_scalar.
// MIT notice and pinned component provenance are retained in inst/COPYRIGHTS.
#include "bed_reader.h"
#include <filesystem>
#include <limits>
#include <stdexcept>

namespace gsim::native {
BedReader::BedReader(const std::string& path, std::uint64_t samples,
                     std::uint64_t markers)
    : path_(path), samples_(samples), markers_(markers),
      stride_(samples / 4u + (samples % 4u != 0u)) {
    const auto fail = [&](const std::string& detail) {
        throw std::runtime_error("BED '" + path_ + "': " + detail);
    };
    if (path.empty() || !samples || !markers) fail("empty path or dimensions");
    if (markers > (std::numeric_limits<std::uint64_t>::max() - 3u) / stride_)
        fail("dimensions overflow file size");
    const auto expected = 3u + markers * stride_;
    if (expected > static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max()) ||
        stride_ > static_cast<std::uint64_t>(std::numeric_limits<std::streamsize>::max()) ||
        stride_ > std::numeric_limits<std::size_t>::max())
        fail("dimensions exceed platform I/O limits");
    const auto file = std::filesystem::u8path(path);
    std::error_code ec;
    const auto actual = std::filesystem::file_size(file, ec);
    if (ec) fail("cannot inspect file: " + ec.message());
    if (actual != expected) fail("file size mismatch (expected " +
        std::to_string(expected) + ", found " + std::to_string(actual) + ")");
    stream_.open(file, std::ios::binary);
    if (!stream_) fail("cannot open file");
    std::uint8_t header[3]{};
    stream_.read(reinterpret_cast<char*>(header), 3);
    if (!stream_ || header[0] != 0x6c || header[1] != 0x1b || header[2] != 1)
        fail("invalid header or non SNP-major mode");
    record_.resize(static_cast<std::size_t>(stride_));
}
void BedReader::read_record(std::uint64_t marker) {
    if (marker >= markers_) throw std::runtime_error("BED '" + path_ + "': marker out of range");
    const auto offset = 3u + marker * stride_; // constructor checked the complete extent
    stream_.clear();
    stream_.seekg(static_cast<std::streamoff>(offset), std::ios::beg);
    stream_.read(reinterpret_cast<char*>(record_.data()), static_cast<std::streamsize>(stride_));
    if (!stream_ || stream_.gcount() != static_cast<std::streamsize>(stride_))
        throw std::runtime_error("BED '" + path_ + "': incomplete record at physical marker " +
                                 std::to_string(marker + 1u));
}
int BedReader::dosage(std::uint64_t sample) const noexcept {
    static const int values[4] = {2, -1, 1, 0}; // BIM A1; -1 is internal missing
    return values[(record_[sample / 4u] >> (2u * (sample % 4u))) & 3u];
}
void BedReader::add_variant(const int* rows, std::uint64_t n, double replacement,
                           double center, double scale, const double* effects,
                           std::uint64_t marker, std::uint64_t markers,
                           std::uint64_t traits, double* output) const {
    // gbits xmat_scalar loop, with caller-frozen gsim sample SD/mean semantics.
    for (std::uint64_t sample = 0; sample < n; ++sample) {
        const int value = dosage(static_cast<std::uint64_t>(rows[sample] - 1));
        const double raw = value < 0 ? replacement : static_cast<double>(value);
        const double transformed = (raw - center) / scale;
        for (std::uint64_t trait = 0; trait < traits; ++trait)
            output[trait * n + sample] += transformed * effects[trait * markers + marker];
    }
}
}
