#ifndef FUMI_TOOLS_CAST_HPP
#define FUMI_TOOLS_CAST_HPP

#include <type_traits>

namespace fumi_tools {

/**
 * Converts the passed type to its corresponding signed type, if necessary.
 */
template <class T>
typename std::make_signed<T>::type as_signed(T t) {
  return static_cast<typename std::make_signed<T>::type>(t);
}

/**
 * Converts the passed type to its corresponding unsigned type, if necessary.
 */
template <class T>
typename std::make_unsigned<T>::type as_unsigned(T t) {
  return static_cast<typename std::make_unsigned<T>::type>(t);
}

}

#endif  // FUMI_TOOLS_CAST_HPP
