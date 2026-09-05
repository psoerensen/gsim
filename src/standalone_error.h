#ifndef GSIM_STANDALONE_ERROR_HPP
#define GSIM_STANDALONE_ERROR_HPP

#include <stdexcept>
#include <string>

namespace gsim::native {

enum Status : int {
    status_success = 0,
    status_invalid_argument = 1,
    status_io_error = 2,
    status_bed_format_error = 3,
    status_out_of_range = 4,
    status_internal_error = 6
};

class Error final : public std::runtime_error {
public:
    Error(Status status, const std::string& message)
        : std::runtime_error(message), status_(status) {}
    Status status() const noexcept { return status_; }
private:
    Status status_;
};

[[noreturn]] inline void invalid_argument(const std::string& message) {
    throw Error(status_invalid_argument, message);
}

} // namespace gsim::native

// Local aliases keep the imported MIT-licensed implementation readable while
// ensuring the implementation itself lives only in gsim::native.
#define GBITS_STATUS_INVALID_ARGUMENT ::gsim::native::status_invalid_argument
#define GBITS_STATUS_IO_ERROR ::gsim::native::status_io_error
#define GBITS_STATUS_BED_FORMAT_ERROR ::gsim::native::status_bed_format_error
#define GBITS_STATUS_OUT_OF_RANGE ::gsim::native::status_out_of_range

#endif
