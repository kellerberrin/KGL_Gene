//
// Created by kellerberrin on 12/11/17.
//

#ifndef KGL_SEQUENCE_VIRTUAL_H
#define KGL_SEQUENCE_VIRTUAL_H

#include <string>
#include <set>
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
  explicit AlphabetSequence(const AlphabetSequence&) = default;
  AlphabetSequence() = delete;
  ~AlphabetSequence() override = default;

  auto operator[] (ContigOffset_t& offset) const { return alphabet_string_[offset]; }
  auto at(ContigOffset_t& offset) const { return alphabet_string_[offset]; }

  ContigSize_t length() const { return alphabet_string_.length(); }
  std::string getSequenceAsString() const override { return alphabet_string_.str(); }

  const std::string compareSequences(const std::shared_ptr<const AlphabetSequence>& sequence_ptr) const {

    return compareSequences(*sequence_ptr, *this);

  }

protected:

  AlphabetString<Alphabet> alphabet_string_;


  static std::string emphasizeDifferences(const AlphabetSequence& alphabet_sequence,
                                          const std::vector<ContigOffset_t>& emphasize_offsets);

  static std::string compareSequences(const AlphabetSequence& compare_sequence,
                                      const AlphabetSequence& reference_sequence);

};


template<typename Alphabet>
std::string AlphabetSequence<Alphabet>::compareSequences(const AlphabetSequence<Alphabet>& compare_sequence,
                                                         const AlphabetSequence<Alphabet>& reference_sequence) {

  ContigSize_t length_compare = compare_sequence.length();
  ContigSize_t length_reference = reference_sequence.length();

  ContigSize_t common_length = length_reference <= length_compare ? length_reference : length_compare;

  std::vector<ContigOffset_t> differences;
  for (ContigOffset_t index = 0; index < common_length; ++index) {


    if (reference_sequence[index] != compare_sequence[index]) {

      differences.push_back(index);

    }

  }

  ExecEnv::log().info("compareSequences(), found: {} sequence differences, common length: {}",
                      differences.size(), common_length);

  std::string comparison_string;
  comparison_string = AlphabetSequence<Alphabet>::emphasizeDifferences(compare_sequence, differences);
  if (common_length < length_compare) {

    comparison_string += "  comparison>";
    for (ContigOffset_t index = common_length; index < length_compare; ++index) {

      comparison_string += Alphabet::convertToChar(compare_sequence[index]);

    }


  } else if (common_length < length_reference) {

    comparison_string += "  reference>";
    for (ContigOffset_t index = common_length; index < length_reference; ++index) {

      comparison_string += Alphabet::convertToChar(reference_sequence[index]);

    }

  }

  return comparison_string;

}

template<typename Alphabet>
std::string AlphabetSequence<Alphabet>::emphasizeDifferences(const AlphabetSequence<Alphabet>& alphabet_sequence,
                                                             const std::vector<ContigOffset_t>& emphasize_offsets) {

  std::string sequence_string = alphabet_sequence.getSequenceAsString();

  if (emphasize_offsets.empty()) return sequence_string;

  // Order the offsets. A < B < C
  std::set<ContigOffset_t> ordered_offsets;
  for (auto offset : emphasize_offsets) {

    ordered_offsets.insert(offset);

  }

  std::string emph_sequence_string;
  size_t index = 0;
  for (auto offset : ordered_offsets) {

    if (offset >= sequence_string.length()) {

      ExecEnv::log().error("emphasizeDifferences() emphasize offset: {} >= sequence string length: {}",
                           offset, emph_sequence_string.length());
      break;
    }

    while (index < offset) {

      emph_sequence_string += sequence_string[index];
      index++;

    }

    // Add spaces to emphasize.
    if (index == offset) {

      emph_sequence_string += ' ';
      emph_sequence_string += std::tolower(sequence_string[index]);
      emph_sequence_string += ' ';
      index++;

    }

  }

  // Add in the rest of the sequence
  while (index < sequence_string.length()) {

    emph_sequence_string += sequence_string[index];
    index++;

  }

  return emph_sequence_string;

}




}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_VIRTUAL_H
