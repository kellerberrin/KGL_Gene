//
// Created by kellerberrin on 12/11/17.
//

#ifndef KGL_SEQUENCE_VIRTUAL_H
#define KGL_SEQUENCE_VIRTUAL_H

#include <string>


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A virtual Sequence class to return the sequences of the derived DNA5Sequence and AminoSequence classes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class AlphabetSequence {

public:

  AlphabetSequence() = default;
  virtual ~AlphabetSequence() = default;

  virtual std::string getSequenceAsString() const = 0;

private:


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_VIRTUAL_H
