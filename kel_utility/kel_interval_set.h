//
// Created by kellerberrin on 17/06/23.
//

#ifndef KEL_INTERVAL_SET_H
#define KEL_INTERVAL_SET_H

#include "kel_interval_unsigned.h"


namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::set.
class IntervalSetLower : public std::set<OpenRightUnsigned, CompareUnsignedIntervalLower> {

public:

  IntervalSetLower() = default;
  ~IntervalSetLower() = default;

  [[nodiscard]] bool containsInterval(const OpenRightUnsigned &interval) const;
  // Returns a vector (possibly empty) of set intervals that intersect the argument interval.
  [[nodiscard]] std::vector<OpenRightUnsigned> findIntersectsInterval(const OpenRightUnsigned &interval) const;
  // Returns a bool the interval and set intersect.
  [[nodiscard]] bool intersectsInterval(const OpenRightUnsigned &interval) const { return not findIntersectsInterval(interval).empty(); }
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
class IntervalSetUpper : public std::set<OpenRightUnsigned, CompareUnsignedIntervalUpper> {

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
class IntervalMultiSetLower : public std::multiset<OpenRightUnsigned, CompareUnsignedIntervalLower> {

public:

  IntervalMultiSetLower() = default;
  ~IntervalMultiSetLower() = default;

};


} // namespace



#endif //KEL_INTERVAL_SET_H
