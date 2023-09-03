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
class IntervalSetLower : public std::set<OpenRightInterval, CompareIntervalLower> {

public:

  IntervalSetLower() = default;
  ~IntervalSetLower() = default;

  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const;
  // Returns a vector (possibly empty) of set intervals that intersect the argument interval.
  [[nodiscard]] std::vector<OpenRightInterval> findIntersectsInterval(const OpenRightInterval &interval) const;
  // Returns a bool the interval and set intersect.
  [[nodiscard]] bool intersectsInterval(const OpenRightInterval &interval) const { return not findIntersectsInterval(interval).empty(); }
  // Returns the interval unions between this set and the argument set.
  // Intervals are modified/extended as necessary for a disjoint interval minimal union.
  [[nodiscard]] IntervalSetLower intervalSetUnion(const IntervalSetLower& interval_set) const;
  // Set with simplified intervals which are modified/extended as necessary for disjoint representation.
  [[nodiscard]] IntervalSetLower simplifyDisjoint() const { return intervalSetUnion(IntervalSetLower()); }


};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::set.
class IntervalSetUpper : public std::set<OpenRightInterval, CompareIntervalUpper> {

public:

  IntervalSetUpper() = default;
  ~IntervalSetUpper() = default;


};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::multiset.
class IntervalMultiSetLower : public std::multiset<OpenRightInterval, CompareIntervalLower> {

public:

  IntervalMultiSetLower() = default;
  ~IntervalMultiSetLower() = default;

};


} // namespace



#endif //KEL_INTERVAL_SET_H
