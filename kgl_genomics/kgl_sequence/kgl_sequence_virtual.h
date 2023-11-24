//
// Created by kellerberrin on 12/11/17.
//

#ifndef KGL_SEQUENCE_VIRTUAL_H
#define KGL_SEQUENCE_VIRTUAL_H


#include "kgl_alphabet_string.h"
#include "kel_interval_unsigned.h"

#include <string>
#include <set>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A virtual Sequence class returns char sequences of the derived DNA5Sequence and AminoSequence classes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VirtualSequence {

public:

  VirtualSequence() = default;
  virtual ~VirtualSequence() = default;

  [[nodiscard]] virtual std::string_view getStringView() const = 0;

};




}   // end namespace

#endif //KGL_SEQUENCE_VIRTUAL_H
