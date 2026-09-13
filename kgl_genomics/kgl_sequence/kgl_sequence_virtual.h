//
// Created by kellerberrin on 12/11/17.
//

#ifndef KGL_SEQUENCE_VIRTUAL_H
#define KGL_SEQUENCE_VIRTUAL_H


#include "kel_interval_unsigned.h"
#include "kel_search.h"

#include <string_view>
#include <vector>


namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A virtual Sequence class returns char sequences of the derived DNA5Sequence and AminoSequence classes.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VirtualSequence {

public:

  VirtualSequence() = default;
  virtual ~VirtualSequence() = default;

  /// ReturnType the sequence as a std::string_view. Assumes that sequence alphabets are 1 byte (char).
  [[nodiscard]] virtual std::string_view getStringView() const = 0;
  /// Searches the underlying sequence using a regular expression. Returns a vector (possibly empty) of intervals.
  [[nodiscard]] std::vector<OpenRightUnsigned> regexSearch(const std::string_view& regex) const {

    return Search::searchView(regex, getStringView());

  }

};


}   // end namespace

#endif //KGL_SEQUENCE_VIRTUAL_H