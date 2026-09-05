#ifndef NATIVE_METADATA_STATUS_HPP
#define NATIVE_METADATA_STATUS_HPP

#include <string>
#include <utility>

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

class Status {
 public:
  Status() = default;
  Status(StatusCode code, std::string message)
      : code_(code), message_(std::move(message)) {}

  static Status ok() { return {}; }

  [[nodiscard]] bool is_ok() const noexcept {
    return code_ == StatusCode::success;
  }
  explicit operator bool() const noexcept { return is_ok(); }
  [[nodiscard]] StatusCode code() const noexcept { return code_; }
  [[nodiscard]] const std::string& message() const noexcept { return message_; }

 private:
  StatusCode code_{StatusCode::success};
  std::string message_{};
};

}  // namespace gsim::native::metadata

#endif

