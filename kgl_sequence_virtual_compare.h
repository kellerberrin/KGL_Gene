//
// Created by kellerberrin on 22/12/17.
//

#ifndef KGL_SEQUENCE_VIRTUAL_COMPARE_H
#define KGL_SEQUENCE_VIRTUAL_COMPARE_H


#include "kgl_sequence_virtual.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



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

    comparison_string += " +(";

    for (ContigOffset_t index = common_length; index < length_compare; ++index) {

      comparison_string += std::tolower(Alphabet::convertToChar(compare_sequence[index]));

    }

    comparison_string += ")";

  } else if (common_length < length_reference) {

    comparison_string += " +(";
    for (ContigOffset_t index = common_length; index < length_reference; ++index) {

      comparison_string += Alphabet::convertToChar(reference_sequence[index]);

    }

    comparison_string += ")";

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


#endif //KGL_SEQUENCE_VIRTUAL_COMPARE_H
