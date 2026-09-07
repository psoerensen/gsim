#ifndef GSIM_BED_READER_H
#define GSIM_BED_READER_H

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace gsim::native {
// Minimal scalar BED reader adapted from gbits (see inst/COPYRIGHTS).
// One packed physical record; no resident genotype panel or external ABI.
class BedReader final {
public:
    BedReader(const std::string& path, std::uint64_t samples, std::uint64_t markers);
    void read_record(std::uint64_t marker);
    int dosage(std::uint64_t sample) const noexcept;
    void add_variant(const int* rows, std::uint64_t n, double replacement,
                     double center, double scale, const double* effects,
                     std::uint64_t marker, std::uint64_t markers,
                     std::uint64_t traits, double* output) const;
private:
    std::string path_;
    std::ifstream stream_;
    std::uint64_t samples_, markers_, stride_;
    std::vector<std::uint8_t> record_;
};
}
#endif
