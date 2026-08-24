//
// Created by kellerberrin on 11/12/23.
//

#ifndef KEL_SEARCH_H
#define KEL_SEARCH_H

#include <regex>
#include "kel_interval_unsigned.h"

namespace kellerberrin {   //  organization level namespace


// Prefer namespace to object.
// Note: std::regex objects are not guaranteed to be thread-safe. The
// const std::regex& overload should only be called from multiple threads
// if the caller provides external synchronization for the regex object.
namespace Search {

  /// Returns a vector of intervals where the regex pattern was found in the view.
  [[nodiscard]] std::vector<OpenRightUnsigned> searchView( const std::regex& search_spec, std::string_view sequence_view);
  /// Convenience routine with regular search expression as text. Use for single ad-hoc text searches.
  [[nodiscard]] std::vector<OpenRightUnsigned> searchView( std::string_view search_spec, std::string_view sequence_view);


};


} // Name space

#endif //KEL_SEARCH_H
