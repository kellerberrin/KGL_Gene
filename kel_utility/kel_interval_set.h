//
// Created by kellerberrin on 17/06/23.
//

#ifndef KEL_INTERVAL_SET_H
#define KEL_INTERVAL_SET_H

#include "kel_interval.h"


namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::set.
class IntervalSet : public std::set<OpenRightInterval, CompareInterval> {

public:

  IntervalSet() = default;
  ~IntervalSet() = default;

  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const;
  [[nodiscard]] bool containsOffset(size_t offset) const { return contains(OpenRightInterval(offset, offset + 1)); }
  // Returns the interval intersection between this set and the argument set returned as an interval set (possibly empty).
  [[nodiscard]] IntervalSet intervalSetIntersection(const IntervalSet& interval_set) const;
  // Returns a vector (possibly empty) of set intervals that intersect the argument interval.
  [[nodiscard]] std::vector<OpenRightInterval> findIntersectsInterval(const OpenRightInterval &interval) const;
  // Returns a bool the interval and set intersect.
  [[nodiscard]] bool intersectsInterval(const OpenRightInterval &interval) const { return not findIntersectsInterval(interval).empty(); }
  // Returns the interval unions between this set and the argument set.
  // Intervals are modified/extended as necessary for a disjoint interval minimal union.
  [[nodiscard]] IntervalSet intervalSetUnion(const IntervalSet& interval_set) const;
  // Set with simplified intervals which are modified/extended as necessary for disjoint representation.
  [[nodiscard]] IntervalSet simplifyDisjoint() const { return intervalSetUnion(IntervalSet()); }


};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::multiset.
class IntervalMultiSet : public std::multiset<OpenRightInterval, CompareInterval> {

public:

  IntervalMultiSet() = default;
  ~IntervalMultiSet() = default;

};


} // namespace



#endif //KEL_INTERVAL_SET_H
