//
// Created by kellerberrin on 15/06/23.
//

#ifndef KEL_INTERVAL_UNSIGNED_H
#define KEL_INTERVAL_UNSIGNED_H


#include <set>
#include <map>
#include <vector>
#include <string>

#include "kel_interval_type.h"

namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Define the unsigned simple right open interval [uint64_t lower_, uint64_t upper_)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

using OpenRightUnsigned = OpenRightInterval<uint64_t, int64_t>;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::set.
class IntervalSetLower : public std::set<OpenRightUnsigned> {

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

using IntervalMultiSetLower = std::multiset<OpenRightUnsigned>;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



struct CompareUnsignedIntervalUpper {

  bool operator()(const OpenRightUnsigned &lhs, const OpenRightUnsigned &rhs) const { return lhs.upper() < rhs.upper(); }

};


// Interval adapter for std::multimap.
template<typename ValueType>
using IntervalUpperMultiMapType = std::multimap<OpenRightUnsigned, ValueType, CompareUnsignedIntervalUpper>;


} // namespace



#endif //KEL_INTERVAL_UNSIGNED_H
