#ifndef GSIM_STANDALONE_METADATA_ERROR_HPP
#define GSIM_STANDALONE_METADATA_ERROR_HPP

#include <stdexcept>
#include <string>

#include "standalone_metadata_status.h"

namespace gsim::native::metadata {

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
