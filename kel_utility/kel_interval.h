//
// Created by kellerberrin on 15/06/23.
//

#ifndef KEL_INTERVAL_H
#define KEL_INTERVAL_H


#include <set>
#include <map>
#include <vector>
#include <string>


namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Defines a simple right open interval [lower_, upper_)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


class OpenRightUnsigned {

public:

  OpenRightUnsigned(size_t lower, size_t upper) { resize(lower, upper); }
  ~OpenRightUnsigned() = default;
  OpenRightUnsigned(const OpenRightUnsigned &copy) = default;

  OpenRightUnsigned &operator=(const OpenRightUnsigned &copy) = default;

  void resize(size_t lower, size_t upper); // Change the interval.
  // Return the shifted interval without changing it's size. Interval is not translated if lower() < 0.
  [[nodiscard]] OpenRightUnsigned translate(int64_t shift) const;
  // Return the zero-translated interval so that lower() == 0.
  [[nodiscard]] OpenRightUnsigned translateZero() const;

  [[nodiscard]] size_t lower() const { return lower_; }
  [[nodiscard]] size_t upper() const { return upper_; }
  [[nodiscard]] size_t size() const { return upper_ - lower_; }

  // Returns the intersection interval or the empty [0, 0) interval indicating no intersection.
  [[nodiscard]] OpenRightUnsigned intersection(const OpenRightUnsigned &interval) const;
  // Merge intersecting or adjacent intervals. If the argument intervals are disjoint and not adjacent
  // then the empty [0, 0) intervals is returned. Note that the merging of the empty intervals will also produce an empty interval.
  [[nodiscard]] OpenRightUnsigned merge(const OpenRightUnsigned &interval) const;

  [[nodiscard]] bool empty() const { return size() == 0; }
  [[nodiscard]] bool containsOffset(size_t offset) const { return offset >= lower_ and offset < upper_; }
  [[nodiscard]] bool containsInterval(const OpenRightUnsigned &interval) const;
  // Note that empty [0, 0) intervals adjoin each other, other empty intervals such as [k, k) and [l, l) do not adjoin if k != l.
  [[nodiscard]] bool adjacent(const OpenRightUnsigned &interval) const { return lower_ == interval.upper_ or interval.lower_ == upper_; }
  [[nodiscard]] bool intersects(const OpenRightUnsigned &interval) const { return not disjoint(interval); }
  [[nodiscard]] bool disjoint(const OpenRightUnsigned &interval) const { return intersection(interval).empty(); }

  // Convenience routine to convert an interval to a string.
  [[nodiscard]] std::string toString() const { return "[ " + std::to_string(lower_) + ", " + std::to_string(upper_) + ")"; }

private:

  size_t upper_{0};
  size_t lower_{0};

};

[[nodiscard]] bool operator==(const OpenRightUnsigned& lhs, const OpenRightUnsigned &rhs);
[[nodiscard]] bool operator<(const OpenRightUnsigned &lhs, const OpenRightUnsigned &rhs);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Comparison operators used to order intervals within indexed containers std::set or std::map.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct CompareIntervalLower {

  bool operator()(const OpenRightUnsigned &lhs, const OpenRightUnsigned &rhs) const { return lhs.lower() < rhs.lower(); }

};


struct CompareIntervalUpper {

  bool operator()(const OpenRightUnsigned &lhs, const OpenRightUnsigned &rhs) const { return lhs.upper() < rhs.upper(); }

};




} // namespace



#endif //KEL_INTERVAL_H
