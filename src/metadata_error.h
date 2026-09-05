#ifndef GSIM_METADATA_ERROR_H
#define GSIM_METADATA_ERROR_H

#include <stdexcept>
#include <string>

namespace gsim::native::metadata {

enum class StatusCode {
  success = 0,
  invalid_argument,
  invalid_extent,
  dimension_mismatch,
  non_finite_value,
  invalid_sparse_structure,
  internal_error
};

class Error : public std::runtime_error {
 public:
  Error(StatusCode code, const std::string& message)
      : std::runtime_error(message), code_(code) {}

  [[nodiscard]] StatusCode code() const noexcept { return code_; }

 private:
  StatusCode code_;
};

}  // namespace gsim::native::metadata

#endif
