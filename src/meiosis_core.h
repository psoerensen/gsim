#pragma once
#include <cstdint>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
namespace meiosis {
inline std::uint64_t mix64(std::uint64_t z) noexcept {
    z = (z ^ (z >> 30u)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27u)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31u);
}

class SplitMix64 {
public:
    explicit SplitMix64(std::uint64_t seed) noexcept : state_(seed) {}

    std::uint64_t next() noexcept {
        state_ += UINT64_C(0x9e3779b97f4a7c15);
        return mix64(state_);
    }

    double open01() noexcept {
        const std::uint64_t leading53 = next() >> 11u;
        return (static_cast<double>(leading53) + 0.5) *
               (1.0 / 9007199254740992.0);
    }

private:
    std::uint64_t state_;
};

inline std::uint64_t stream_seed_raw(std::uint64_t seed, std::uint64_t animal, std::uint64_t chromosome,
                          int side) {
    std::uint64_t key = mix64(seed ^ UINT64_C(0x243f6a8885a308d3));
    key ^= mix64(animal ^
                 UINT64_C(0x13198a2e03707344));
    key ^= mix64(chromosome ^
                 UINT64_C(0xa4093822299f31d0));
    key ^= side == 1 ? UINT64_C(0x082efa98ec4e6c89)
                     : UINT64_C(0x452821e638d01377);
    return mix64(key);
}

struct Events {
    std::vector<double> crossovers;
    int starting_haplotype;
};

inline Events draw_events_raw(const double* positions, std::size_t count, std::uint64_t seed, std::uint64_t animal,
                   std::uint64_t chromosome, int side) {
    SplitMix64 rng(stream_seed_raw(seed, animal, chromosome, side));
    const double first = positions[0];
    const double last = positions[count - 1];
    const double length = last - first;
    std::size_t crossover_count = 0u;
    if (length > 0.0) {
        double arrival = 0.0;
        for (;;) {
            arrival += -std::log(rng.open01());
            if (arrival > length) break;
            if (crossover_count == std::numeric_limits<std::size_t>::max()) {
                throw std::runtime_error("crossover count exceeds the supported range");
            }
            ++crossover_count;
        }
    }
    std::vector<double> crossovers(crossover_count);
    for (double& crossover : crossovers) {
        crossover = first + length * rng.open01();
    }
    std::sort(crossovers.begin(), crossovers.end());
    const int starting_haplotype = rng.open01() < 0.5 ? 1 : 2;
    return {std::move(crossovers), starting_haplotype};
}

}
