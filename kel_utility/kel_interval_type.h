//
// Created by kellerberrin on 21/09/23.
//

#ifndef KEL_INTERVAL_TYPE_H
#define KEL_INTERVAL_TYPE_H


#include <set>
#include <map>
#include <vector>
#include <string>

#include "kel_exec_env.h"


namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Define the concepts for the underlying numeric types used in the OpenRightInterval template class.
// This includes floats. However, floats may require work on the comparison operators.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
concept IntervalNumericType = std::integral<T> || std::floating_point<T>;

template < class T >
concept SignedIntervalType = (std::integral<T> && std::is_signed_v<T>) || std::floating_point<T>;

template < class T >
concept UnsignedIntervalType = std::integral<T> && !std::signed_integral<T>;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Defines a simple right open interval [lower_, upper_)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <IntervalNumericType IntervalValue, SignedIntervalType SignedInterval>
class OpenRightInterval {

public:

  constexpr OpenRightInterval(IntervalValue lower, IntervalValue upper) { resize(lower, upper); }
  constexpr ~OpenRightInterval() = default;
  constexpr OpenRightInterval(const OpenRightInterval &copy) = default;

  constexpr OpenRightInterval &operator=(const OpenRightInterval &copy) = default;

  constexpr void resize(IntervalValue lower, IntervalValue upper)  {

    if (upper < lower) {

      ExecEnv::log().warn("OpenRightInterval::resize; Incorrect Initialization, Upper Offset: {} < Lower Offset: {}", upper, lower);
      std::swap(lower, upper);

    }

    lower_ = lower;
    upper_ = upper;

  }


  // Shift the interval without changing it's size.
  [[nodiscard]] constexpr OpenRightInterval translate(SignedInterval shift) const {

    OpenRightInterval translated(*this);

    if constexpr (UnsignedIntervalType<IntervalValue>) {

      if ((static_cast<SignedInterval>(lower()) + shift) < 0) {

        ExecEnv::log().warn("OpenRightInterval::translate; translate shift: {} results in negative values for interval: {}", shift, toString());
        return {lower(), upper()};

      }

    }

    IntervalValue translate_lower = static_cast<SignedInterval>(lower()) + shift;
    IntervalValue translate_upper = static_cast<SignedInterval>(upper()) + shift;

    return { translate_lower , translate_upper };

  }

  // Return the zero-translated interval so that lower() == 0.
  [[nodiscard]] constexpr OpenRightInterval translateZero() const {

    SignedInterval shift =  -1 * static_cast<SignedInterval>(lower());
    return translate(shift);

  }

  [[nodiscard]] constexpr IntervalValue lower() const { return lower_; }

  [[nodiscard]] constexpr IntervalValue upper() const { return upper_; }

  [[nodiscard]] constexpr size_t size() const { return upper_ - lower_; }

  // Returns the intersection interval or the empty [0, 0) interval indicating no intersection.
  [[nodiscard]] constexpr OpenRightInterval intersection(const OpenRightInterval &interval) const {

    if (lower_ >= interval.upper_ or interval.lower_ >= upper_) {

      return {0, 0};

    }

    return { std::max<IntervalValue>(lower_, interval.lower_), std::min<IntervalValue>(upper_, interval.upper_ ) };

  }

  // Merge intersecting or adjacent intervals. If the argument intervals are disjoint and not adjacent
  // then the empty [0, 0) interval is returned.
  // Note that the merging of the empty intervals will also produce an empty interval.
  [[nodiscard]] constexpr OpenRightInterval merge(const OpenRightInterval &interval) const {

    if (intersects(interval) or adjacent(interval)) {

      return { std::min<IntervalValue>(lower_, interval.lower_), std::max<IntervalValue>(upper_, interval.upper_ ) };

    }

    return {0, 0};

  }

  [[nodiscard]] constexpr bool empty() const { return size() == 0; }

  [[nodiscard]] constexpr bool containsOffset(size_t offset) const { return offset >= lower_ and offset < upper_; }

  [[nodiscard]] constexpr bool containsInterval(const OpenRightInterval &interval) const { return intersection(interval) == interval; }

  // Note that empty [0, 0) intervals adjoin each other, other empty intervals such as [k, k) and [l, l) do not adjoin if k != l.
  [[nodiscard]] constexpr bool adjacent(const OpenRightInterval &interval) const { return lower_ == interval.upper_ or interval.lower_ == upper_; }

  [[nodiscard]] constexpr bool intersects(const OpenRightInterval &interval) const { return not disjoint(interval); }

  [[nodiscard]] constexpr bool disjoint(const OpenRightInterval &interval) const { return intersection(interval).empty(); }

  // Convenience routine to convert an interval to a string.
  [[nodiscard]] constexpr std::string toString() const { return "[ " + std::to_string(lower_) + ", " + std::to_string(upper_) + ")"; }

  // Define an ordering using the spaceship operator. Intervals are, by default, ordered by their lower value.
  constexpr auto operator<=>(const OpenRightInterval &rhs) const { return std::tie(lower_, upper_) <=> std::tie(rhs.lower_, rhs.upper_); }
  constexpr bool operator==(const OpenRightInterval &rhs) const { return lower() == rhs.lower() and upper() == rhs.upper(); }

private:

  IntervalValue upper_{0};
  IntervalValue lower_{0};

};


} // namespace


#endif //KEL_INTERVAL_TYPE_H
