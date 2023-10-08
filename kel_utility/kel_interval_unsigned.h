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
// Define the unsigned simple right open interval [size_t lower_, size_t upper_)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

using OpenRightUnsigned = OpenRightInterval<size_t, int64_t>;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Comparison operators used to order intervals within indexed containers std::set or std::map.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct CompareUnsignedIntervalLower {

  bool operator()(const OpenRightUnsigned &lhs, const OpenRightUnsigned &rhs) const { return lhs.lower() < rhs.lower(); }

};


struct CompareUnsignedIntervalUpper {

  bool operator()(const OpenRightUnsigned &lhs, const OpenRightUnsigned &rhs) const { return lhs.upper() < rhs.upper(); }

};




} // namespace



#endif //KEL_INTERVAL_UNSIGNED_H
