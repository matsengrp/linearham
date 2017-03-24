#ifndef BRENT_FUNC_
#define BRENT_FUNC_

#include <brent.hpp>

/// @file brent_func.hpp
/// @brief Header for the Brent optimization functor template.

namespace brent {

/// @brief An alias template for (input/output: double) member function
/// pointers.
template <class T>
using DDMemFnPtr = double (T::*)(double);


/// @brief A class template for Brent optimization functors based on
/// (input/output: double) member functions.
template <class T>
class member_func_wrapper : public func_base {
 private:
  T* obj_;
  DDMemFnPtr<T> func_;

 public:
  member_func_wrapper(){};
  member_func_wrapper(T* obj, DDMemFnPtr<T> func) : obj_(obj), func_(func){};
  double operator()(double x) override { return (obj_->*func_)(x); };
};
}

#endif  // BRENT_FUNC_
