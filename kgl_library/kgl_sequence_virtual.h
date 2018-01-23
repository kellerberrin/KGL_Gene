//
// Created by kellerberrin on 12/11/17.
//

#ifndef KGL_SEQUENCE_VIRTUAL_H
#define KGL_SEQUENCE_VIRTUAL_H

#include <string>
#include <set>
#include "kgl_alphabet_string.h"
#include "kgl_sequence_manip.h"



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
  AlphabetSequence() = default;
  ~AlphabetSequence() override = default;

  auto operator[] (ContigOffset_t offset) const { return alphabet_string_[offset]; }
  auto at(ContigOffset_t offset) const { return alphabet_string_[offset]; }

  ContigSize_t length() const { return alphabet_string_.length(); }
  std::string getSequenceAsString() const override { return alphabet_string_.str(); }
  const AlphabetString<Alphabet>& getAlphabetString() const { return alphabet_string_; }

protected:

  AlphabetString<Alphabet> alphabet_string_;

//  const std::string compareSequences(const AlphabetSequence& compare_seq) const {  return compareSequences(compare_seq, *this); }

  const std::string compareSequencesDNA(const AlphabetSequence& compare, CompareScore_t& score) const;
  const std::string compareSequencesAmino(const AlphabetSequence& compare, CompareScore_t& score) const;
  CompareScore_t scoreMyerHirschberg(const AlphabetSequence& compare) const;
  CompareScore_t scoreLevenshtein(const AlphabetSequence<Alphabet>& compare_sequence) const;
  static const std::string multipleAlign(const std::vector<std::shared_ptr<const AlphabetSequence>>& compare_vec);

  // Letter offset is relative to the begining of the sequence (0 is the first letter).
  bool modifyLetter(ContigOffset_t sequence_offset, typename Alphabet::Alphabet letter);
  // Delete offset is relative to the begining of the sequence (0 is the first letter).
  bool deleteOffset(ContigOffset_t delete_offset, ContigSize_t delete_size);
  // Insert offset is relative to the begining of the sequence (0 is the first letter).
  bool insertOffset(ContigOffset_t insert_offset, const AlphabetSequence& inserted_sequence);
  // Returns bool false if offset and/or size are out of bounds.
  bool getSubsequence(ContigOffset_t substring_offset, ContigSize_t delete_size, std::shared_ptr<AlphabetSequence> sub_sequence) const;

private:


};


template<typename Alphabet>
bool AlphabetSequence<Alphabet>::modifyLetter(ContigOffset_t sequence_offset, typename Alphabet::Alphabet letter) {

  if (sequence_offset >= length()) {

    ExecEnv::log().error("modifyLetter(), sequence offset: {} exceeds sequence size: {}",
                         sequence_offset, length());
    return false;
  }

  alphabet_string_.modifyLetter(sequence_offset, letter);

  return true;

}


template<typename Alphabet>
bool AlphabetSequence<Alphabet>::deleteOffset(ContigOffset_t delete_offset, ContigSize_t delete_size) {

  if ((delete_offset + delete_size) > length()) {

    ExecEnv::log().error("Attempt to delete past the end a sequence string, offset: {}, delete size: {}, sequence size: {}",
                         delete_offset, delete_size, length());
    return false;

  }

  if (not alphabet_string_.erase(delete_offset, delete_size)) {

    ExecEnv::log().error("Problem deleting subsequence from sequence string, offset: {}, delete size: {}, sequence size: {}",
                         delete_offset, delete_size, length());
    return false;

  }

  return true;

}


template<typename Alphabet>
bool AlphabetSequence<Alphabet>::insertOffset(ContigOffset_t insert_offset, const AlphabetSequence& inserted_sequence) {

  if (insert_offset >= length()) {

    ExecEnv::log().error("Attempt to insert past the end a sequence string, insert offset: {}, insert sub-sequence: {}",
                         insert_offset, inserted_sequence.getSequenceAsString());
    return false;

  }

  if (not alphabet_string_.insert(insert_offset, inserted_sequence.alphabet_string_)) {

    ExecEnv::log().error("Problem inserting a sub-sequence , insert offset: {}, insert sub-sequence: {}",
                         insert_offset, inserted_sequence.getSequenceAsString());
    return false;

  }

  return true;

}


template<typename Alphabet>
bool AlphabetSequence<Alphabet>::getSubsequence(ContigOffset_t substring_offset, ContigSize_t substring_size, std::shared_ptr<AlphabetSequence> sub_sequence) const {

  if ((substring_offset + substring_size) > length()) {

    ExecEnv::log().error("Attempt to get substring past the end a sequence string, offset: {}, size: {}, sequence size: {}",
                         substring_offset, substring_size, length());
    return false;

  }

  if (not alphabet_string_.substring(substring_offset, substring_size, sub_sequence->alphabet_string_)) {

    ExecEnv::log().error("Problem generating a subsequence , offset: {}, subsequence size: {}, from sequence size: {}",
                         substring_offset, substring_size, alphabet_string_.length());
    return false;

  }

  return true;

}

template<typename Alphabet>
const std::string AlphabetSequence<Alphabet>::compareSequencesDNA(const AlphabetSequence<Alphabet>& compare_sequence, CompareScore_t& score) const
{

  if (length() == 0 or compare_sequence.length() == 0) {

    ExecEnv::log().error("compareSequences(), Cannot compare empty sequences");
    score = 0;
    return "";

  }

  return SequenceManipulation().compareSequencesDNA(getSequenceAsString(), compare_sequence.getSequenceAsString(), score);

}

template<typename Alphabet>
const std::string AlphabetSequence<Alphabet>::compareSequencesAmino(const AlphabetSequence<Alphabet>& compare_sequence, CompareScore_t& score) const
{

  if (length() == 0 or compare_sequence.length() == 0) {

    ExecEnv::log().error("compareSequences(), Cannot compare empty sequences");
    score = 0;
    return "";

  }

  return SequenceManipulation().compareSequencesAmino(getSequenceAsString(), compare_sequence.getSequenceAsString(), score);

}


template<typename Alphabet>
const std::string AlphabetSequence<Alphabet>::multipleAlign(const std::vector<std::shared_ptr<const AlphabetSequence>>& compare_vec) {

  std::vector<std::string> str_vector;

  for (auto sequence : compare_vec) {

    str_vector.push_back(sequence->getSequenceAsString());

  }

  return SequenceManipulation().compareSequencesMultiple(str_vector);

}


template<typename Alphabet>
CompareScore_t AlphabetSequence<Alphabet>::scoreMyerHirschberg(const AlphabetSequence<Alphabet>& compare_sequence) const {

  if (length() == 0 or compare_sequence.length() == 0) {

    ExecEnv::log().error("compareSequences(), Cannot compare empty sequences");
    return 0;

  }

  return SequenceManipulation().compareMyerHirschberg(getSequenceAsString(), compare_sequence.getSequenceAsString());

}


template<typename Alphabet>
CompareScore_t AlphabetSequence<Alphabet>::scoreLevenshtein(const AlphabetSequence<Alphabet>& compare_sequence) const {

  if (length() == 0 or compare_sequence.length() == 0) {

    ExecEnv::log().error("compareSequences(), Cannot compare empty sequences");
    return 0;

  }

  return SequenceManipulation().Levenshtein(getSequenceAsString(), compare_sequence.getSequenceAsString());

}



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_VIRTUAL_H
