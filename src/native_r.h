#ifndef GSIM_NATIVE_R_H
#define GSIM_NATIVE_R_H

#include <stdexcept>
#include <string>
#include <typeinfo>

#include <R.h>
#include <Rinternals.h>

namespace gsim::native::r {

template <typename T>
SEXP type_tag() {
    return Rf_install(typeid(T).name());
}

template <typename T>
void finalizer(SEXP pointer) noexcept {
    delete static_cast<T*>(R_ExternalPtrAddr(pointer));
    R_ClearExternalPtr(pointer);
}

template <typename T>
T* require(SEXP pointer, const char* description) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == nullptr ||
        R_ExternalPtrTag(pointer) != type_tag<T>()) {
        throw std::runtime_error(std::string(description) +
                                 " are invalid or have been released");
    }
    return static_cast<T*>(R_ExternalPtrAddr(pointer));
}

template <typename T>
SEXP make_owned(T* value) {
    SEXP pointer = PROTECT(R_MakeExternalPtr(value, type_tag<T>(), R_NilValue));
    R_RegisterCFinalizerEx(pointer, finalizer<T>, TRUE);
    UNPROTECT(1);
    return pointer;
}

template <typename T>
void release(SEXP pointer, const char* description) {
    (void)require<T>(pointer, description);
    finalizer<T>(pointer);
}

}  // namespace gsim::native::r

#endif
