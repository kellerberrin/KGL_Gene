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


class OpenRightInterval {

public:

  OpenRightInterval(size_t lower, size_t upper) { resize(lower, upper); }
  ~OpenRightInterval() = default;
  OpenRightInterval(const OpenRightInterval &copy) = default;

  OpenRightInterval &operator=(const OpenRightInterval &copy) = default;

  void resize(size_t lower, size_t upper); // Change the interval.
  void translate(int64_t shift); // Shift the interval without changing it's size.
  // Shift the interval so that lower() == 0. The translation (which is always non-positive) is returned.
  int64_t translateZero();

  [[nodiscard]] size_t lower() const { return lower_; }
  [[nodiscard]] size_t upper() const { return upper_; }
  [[nodiscard]] size_t size() const { return upper_ - lower_; }

  // Returns the intersection interval or the empty [0, 0) interval indicating no intersection.
  [[nodiscard]] OpenRightInterval intersection(const OpenRightInterval &interval) const;
  // Merge intersecting or adjacent intervals. If the argument intervals are disjoint and not adjacent
  // then the empty [0, 0) intervals is returned. Note that the merging of the empty intervals will also produce an empty interval.
  [[nodiscard]] OpenRightInterval merge(const OpenRightInterval &interval) const;

  [[nodiscard]] bool empty() const { return size() == 0; }
  [[nodiscard]] bool containsOffset(size_t offset) const { return offset >= lower_ and offset < upper_; }
  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const;
  // Note that empty [0, 0) intervals adjoin each other, other empty intervals such as [k, k) and [l, l) do not adjoin if k != l.
  [[nodiscard]] bool adjacent(const OpenRightInterval &interval) const { return lower_ == interval.upper_ or interval.lower_ == upper_; }
  [[nodiscard]] bool intersects(const OpenRightInterval &interval) const { return not disjoint(interval); }
  [[nodiscard]] bool disjoint(const OpenRightInterval &interval) const { return intersection(interval).empty(); }

  // Convenience routine to convert an interval to a string.
  [[nodiscard]] std::string toString() const { return "[ " + std::to_string(lower_) + ", " + std::to_string(upper_) + ")"; }

private:

  size_t upper_{0};
  size_t lower_{0};

};

[[nodiscard]] bool operator==(const OpenRightInterval& lhs, const OpenRightInterval &rhs);
[[nodiscard]] bool operator<(const OpenRightInterval &lhs, const OpenRightInterval &rhs);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Comparison operators used to order intervals within indexed containers std::set or std::map.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct CompareIntervalLower {

  bool operator()(const OpenRightInterval &lhs, const OpenRightInterval &rhs) const { return lhs.lower() < rhs.lower(); }

};


struct CompareIntervalUpper {

  bool operator()(const OpenRightInterval &lhs, const OpenRightInterval &rhs) const { return lhs.upper() < rhs.upper(); }

};




} // namespace



#endif //KEL_INTERVAL_H
