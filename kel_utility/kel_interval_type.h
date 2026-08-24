//
// Created by kellerberrin on 21/09/23.
//

#ifndef KEL_INTERVAL_TYPE_H
#define KEL_INTERVAL_TYPE_H


#include <set>
#include <map>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <utility>
#include <limits>

#include "kel_exec_env.h"


namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Define the concepts for the underlying numeric types used in the OpenRightInterval template class.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Ensure the minimum size of the underlying numeric type is 4 bytes.
template <typename T>
concept IntervalNumericType = std::integral<T> && sizeof(T) >= 4;

// Define the signed and unsigned interval types.
template < class T >
concept UnsignedIntervalType = IntervalNumericType<T> && !std::signed_integral<T>;

template < class T >
concept SignedIntervalType = IntervalNumericType<T> && std::is_signed_v<T>;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Defines a simple right open interval [lower_, upper_)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <UnsignedIntervalType IntervalValue>
class OpenRightInterval {

public:

  using SignedInterval = std::make_signed_t<IntervalValue>;

  static_assert(sizeof(IntervalValue) == sizeof(SignedInterval),
                "IntervalValue and its signed counterpart must have the same size");

  constexpr static OpenRightInterval EMPTY_INTERVAL{0, 0};

  constexpr OpenRightInterval(IntervalValue lower, IntervalValue upper) { resize(lower, upper); }
  constexpr ~OpenRightInterval() = default;
  constexpr OpenRightInterval(const OpenRightInterval &copy) = default;

  constexpr OpenRightInterval &operator=(const OpenRightInterval &copy) = default;

  /// Resize the interval to [lower, upper), swapping the bounds if upper < lower.
  constexpr void resize(IntervalValue lower, IntervalValue upper)  {

    if (upper < lower) {

      ExecEnv::log().warn("OpenRightInterval::resize; Incorrect Initialization, Upper Offset: {} < Lower Offset: {}", upper, lower);
      std::swap(lower, upper);

    }

    lower_ = lower;
    upper_ = upper;

  }


  /// Shift the interval without changing its size.
  [[nodiscard]] constexpr OpenRightInterval translate(SignedInterval shift) const {

    OpenRightInterval translated(*this);

    if constexpr (UnsignedIntervalType<IntervalValue>) {

      // Guard against signed overflow when lower() exceeds the signed range.
      if (lower() > static_cast<IntervalValue>(std::numeric_limits<SignedInterval>::max())) {

        ExecEnv::log().warn("OpenRightInterval::translate; lower offset: {} exceeds signed range, cannot translate", lower());
        return {lower(), upper()};

      }

      if ((static_cast<SignedInterval>(lower()) + shift) < 0) {

        ExecEnv::log().warn("OpenRightInterval::translate; translate shift: {} results in negative values for interval: {}", shift, toString());
        return {lower(), upper()};

      }

    }

    IntervalValue translate_lower = static_cast<SignedInterval>(lower()) + shift;
    IntervalValue translate_upper = static_cast<SignedInterval>(upper()) + shift;

    return { translate_lower , translate_upper };

  }

  /// Return the zero-translated interval so that lower() == 0.
  [[nodiscard]] constexpr OpenRightInterval translateZero() const {

    SignedInterval shift =  -1 * static_cast<SignedInterval>(lower());
    return translate(shift);

  }

  /// Return the lower bound of the interval.
  [[nodiscard]] constexpr IntervalValue lower() const { return lower_; }

  /// Return the upper bound of the interval.
  [[nodiscard]] constexpr IntervalValue upper() const { return upper_; }

  /// Return the size (width) of the interval.
  [[nodiscard]] constexpr size_t size() const { return upper_ - lower_; }

  /// Returns the intersection interval or the empty [0, 0) interval indicating no intersection.
  [[nodiscard]] constexpr OpenRightInterval intersection(const OpenRightInterval &interval) const {

    if (lower_ >= interval.upper_ or interval.lower_ >= upper_) {

      return EMPTY_INTERVAL;

    }

    return { std::max<IntervalValue>(lower_, interval.lower_), std::min<IntervalValue>(upper_, interval.upper_ ) };

  }

  /// Merge intersecting or adjacent intervals. If the argument intervals are disjoint and not adjacent
  /// then the empty [0, 0) interval is returned.
  /// Note that the merging of the empty intervals will also produce an empty interval.
  [[nodiscard]] constexpr OpenRightInterval merge(const OpenRightInterval &interval) const {

    if (intersects(interval) or adjacent(interval)) {

      return { std::min<IntervalValue>(lower_, interval.lower_), std::max<IntervalValue>(upper_, interval.upper_ ) };

    }

    return EMPTY_INTERVAL;

  }

  /// Return true if the interval is empty (zero width).
  [[nodiscard]] constexpr bool empty() const { return size() == 0; }

  /// Return true if the offset lies within [lower_, upper_).
  [[nodiscard]] constexpr bool containsOffset(size_t offset) const { return offset >= lower_ and offset < upper_; }

  /// Return true if the argument interval is fully contained within this interval.
  [[nodiscard]] constexpr bool containsInterval(const OpenRightInterval &interval) const { return intersection(interval) == interval; }

  /// Note that empty [0, 0) intervals adjoin each other, other empty intervals such as [k, k) and [l, l) do not adjoin if k != l.
  [[nodiscard]] constexpr bool adjacent(const OpenRightInterval &interval) const { return lower_ == interval.upper_ or interval.lower_ == upper_; }

  /// Return true if the argument interval intersects this interval.
  [[nodiscard]] constexpr bool intersects(const OpenRightInterval &interval) const { return not disjoint(interval); }

  /// Return true if the argument interval is disjoint from this interval.
  [[nodiscard]] constexpr bool disjoint(const OpenRightInterval &interval) const { return intersection(interval).empty(); }

  /// Convenience routine to convert an interval to a string.
  [[nodiscard]] constexpr std::string toString() const { return "[ " + std::to_string(lower_) + ", " + std::to_string(upper_) + ")"; }

  /// Define an ordering using the spaceship operator. Intervals are, by default, ordered by their lower value.
  constexpr auto operator<=>(const OpenRightInterval &rhs) const { return std::tie(lower_, upper_) <=> std::tie(rhs.lower_, rhs.upper_); }
  constexpr bool operator==(const OpenRightInterval &rhs) const { return lower() == rhs.lower() and upper() == rhs.upper(); }


private:

  IntervalValue upper_{0};
  IntervalValue lower_{0};

};


} // namespace


#endif //KEL_INTERVAL_TYPE_H
