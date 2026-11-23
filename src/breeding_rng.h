#pragma once
#include <cstdint>
#include <cmath>
#include <string>
namespace breeding {
inline std::uint64_t mix(std::uint64_t x) {
  x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);return x^(x>>31);
}
inline std::uint64_t hash(const std::string& s) {
  std::uint64_t h=UINT64_C(14695981039346656037);for(unsigned char c:s){h^=c;h*=UINT64_C(1099511628211);}return h;
}
struct Rng {
  using result_type=std::uint64_t;
  std::uint64_t state;
  static constexpr result_type min(){return 0;} static constexpr result_type max(){return UINT64_MAX;}
  result_type operator()(){state+=UINT64_C(0x9e3779b97f4a7c15);return mix(state);}
  double uniform(){return (static_cast<double>((*this)()>>11)+.5)/9007199254740992.;}
  double normal(){return std::sqrt(-2*std::log(uniform()))*std::cos(6.2831853071795864769*uniform());}
};
inline Rng stream(std::uint64_t s,const std::string& id,const std::string& component,const std::string& domain) {
  return {mix(mix(s)^mix(hash(id))^mix(hash(component)^UINT64_C(0x13198a2e03707344))^
    mix(hash(domain)^UINT64_C(0xa4093822299f31d0)))};
}
}
