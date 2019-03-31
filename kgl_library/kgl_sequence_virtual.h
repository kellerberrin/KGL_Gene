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
class  AlphabetSequence : public VirtualSequence {

public:

  explicit AlphabetSequence(AlphabetString<Alphabet>&& sequence) noexcept : alphabet_string_(std::move(sequence)) {}
  AlphabetSequence() = default;
  AlphabetSequence(const AlphabetSequence&) = delete;
  ~AlphabetSequence() override = default;

  AlphabetSequence& operator=(const AlphabetSequence&) = delete;

  auto operator[] (ContigOffset_t offset) const { return alphabet_string_[offset]; }
  auto at(ContigOffset_t offset) const { return alphabet_string_[offset]; }

  ContigSize_t length() const { return alphabet_string_.length(); }
  std::string getSequenceAsString() const override { return alphabet_string_.str(); }
  const AlphabetString<Alphabet>& getAlphabetString() const { return alphabet_string_; }

  // Check for memory corruption.
  bool verifySequence() const { return alphabet_string_.verifyString(); }

  // A hash value of the sequence.
  size_t hashSequence() const { return alphabet_string_.hashString(); }

  // Search for all subsequences.
  std::vector<ContigOffset_t> findAll(const AlphabetSequence& sub_sequence) const { return alphabet_string_.findAll(sub_sequence.alphabet_string_); }


protected:

  AlphabetString<Alphabet> alphabet_string_;

  // Letter offset is relative to the begining of the sequence (0 is the first letter).
  bool modifyLetter(ContigOffset_t sequence_offset, typename Alphabet::Alphabet letter);
  // Delete offset is relative to the begining of the sequence (0 is the first letter).
  bool deleteOffset(ContigOffset_t delete_offset, ContigSize_t delete_size);
  // Insert offset is relative to the begining of the sequence (0 is the first letter).
  bool insertOffset(ContigOffset_t insert_offset, const AlphabetSequence& inserted_sequence);
  // Returns bool false if offset and/or size are out of bounds.
  bool getSubsequence(ContigOffset_t substring_offset, ContigSize_t substring_size, std::shared_ptr<AlphabetSequence> sub_sequence) const;

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

  if (insert_offset < length()) {

    if (not alphabet_string_.insert(insert_offset, inserted_sequence.alphabet_string_)) {

      ExecEnv::log().error("Problem inserting a sub-sequence , insert offset: {}, insert sub-sequence: {}",
                           insert_offset, inserted_sequence.getSequenceAsString());
      return false;

    }

  } else if (insert_offset == length()) {

    if (not alphabet_string_.append(inserted_sequence.alphabet_string_)) {

      ExecEnv::log().error("Problem appending a sub-sequence , insert offset: {}, insert sub-sequence: {}",
                           insert_offset, inserted_sequence.getSequenceAsString());
      return false;

    }

  } else {

    ExecEnv::log().error("Attempt to insert past the end a sequence string length:{}, insert offset: {}, insert sub-sequence: {}",
                         length(), insert_offset, inserted_sequence.getSequenceAsString());
    return false;

  }


  return true;

}


template<typename Alphabet>
bool AlphabetSequence<Alphabet>::getSubsequence(ContigOffset_t substring_offset,
                                                ContigSize_t substring_size,
                                                std::shared_ptr<AlphabetSequence> sub_sequence) const {

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

}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_VIRTUAL_H
