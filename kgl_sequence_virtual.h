//
// Created by kellerberrin on 12/11/17.
//

#ifndef KGL_SEQUENCE_VIRTUAL_H
#define KGL_SEQUENCE_VIRTUAL_H

#include <string>
#include "kgl_alphabet_string.h"



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A virtual Sequence class to return the sequences of the derived DNA5Sequence and AminoSequence classes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VirtualSequence {

public:

  VirtualSequence() = default;
  virtual ~VirtualSequence() = default;

  virtual std::string getSequenceAsString() const = 0;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A template Sequence class to hold the alphabet strings for each sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename Alphabet>
class AlphabetSequence : public VirtualSequence {

public:

  explicit AlphabetSequence(AlphabetString<Alphabet> sequence) : alphabet_string_(std::move(sequence)) {}
  AlphabetSequence() = delete;
  ~AlphabetSequence() override = default;


  typename Alphabet::Alphabet operator[] (ContigOffset_t& offset) const { return alphabet_string_[offset]; }
  typename Alphabet::Alphabet at(ContigOffset_t& offset) const { return alphabet_string_[offset]; }

  ContigSize_t length() const { return alphabet_string_.length(); }

  std::string getSequenceAsString() const override { return alphabet_string_.str(); }



protected:

  AlphabetString<Alphabet> alphabet_string_;


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_VIRTUAL_H
