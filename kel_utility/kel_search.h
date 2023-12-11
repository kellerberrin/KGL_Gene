//
// Created by kellerberrin on 11/12/23.
//

#ifndef KEL_SEARCH_H
#define KEL_SEARCH_H

#include <regex>
#include "kel_interval_unsigned.h"

namespace kellerberrin {   //  organization level namespace


// Object cannot be created, just supplies scope and visibility.
class Search {

public:

  Search() = delete;
  ~Search() = delete;

  // Returns a vector of intervals where the regex pattern was found in the view.
  [[nodiscard]] static std::vector<OpenRightUnsigned> searchView(const std::regex& search_spec, const std::string_view& sequence_view);
  // Convenience routine with regular search expression as text. Use for single ad-hoc text searches.
  [[nodiscard]] static std::vector<OpenRightUnsigned> searchView(const std::string_view& search_spec, const std::string_view& sequence_view);


};


} // Name space

#endif //KEL_SEARCH_H
