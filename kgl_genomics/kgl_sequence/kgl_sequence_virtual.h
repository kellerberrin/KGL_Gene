//
// Created by kellerberrin on 12/11/17.
//

#ifndef KGL_SEQUENCE_VIRTUAL_H
#define KGL_SEQUENCE_VIRTUAL_H

#include <string>
#include <set>
#include "kgl_alphabet_string.h"



namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A virtual Sequence class to return the sequences of the derived DNA5Sequence and AminoSequence classes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VirtualSequence {

public:

  VirtualSequence() = default;
  virtual ~VirtualSequence() = default;

  [[nodiscard]] virtual std::string getSequenceAsString() const = 0;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A template Sequence class to hold the alphabet strings for each sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename Alphabet>
class  AlphabetSequence : public VirtualSequence {

public:

  AlphabetSequence(AlphabetSequence&& sequence) noexcept : alphabet_string_(std::move(sequence.alphabet_string_)) {}
  explicit AlphabetSequence(AlphabetString<Alphabet>&& sequence) noexcept : alphabet_string_(std::move(sequence)) {}
  AlphabetSequence() = default;
  AlphabetSequence(const AlphabetSequence&) = delete; // For Performance reasons, don't allow copy constructors
  ~AlphabetSequence() override = default;

  // For Performance reasons, don't allow naive assignments.
  AlphabetSequence& operator=(const AlphabetSequence&) = delete;
  // Only allow move assignments
  AlphabetSequence& operator=(AlphabetSequence&& moved) noexcept {

    alphabet_string_ = std::move(moved.alphabet_string_);
    return *this;

  }

  [[nodiscard]] auto operator[] (ContigOffset_t offset) const { return alphabet_string_[offset]; }
  [[nodiscard]] auto at(ContigOffset_t offset) const { return alphabet_string_[offset]; }

  [[nodiscard]] ContigSize_t length() const { return alphabet_string_.length(); }
  [[nodiscard]] std::string getSequenceAsString() const override { return alphabet_string_.str(); }
  [[nodiscard]] const AlphabetString<Alphabet>& getAlphabetString() const { return alphabet_string_; }

  // Check for memory corruption.
  [[nodiscard]] bool verifySequence() const { return alphabet_string_.verifyString(); }

  // A hash value of the sequence.
  [[nodiscard]] size_t hashSequence() const { return alphabet_string_.hashString(); }

  // Search for all subsequences.
  [[nodiscard]] std::vector<ContigOffset_t> findAll(const AlphabetSequence& sub_sequence) const { return alphabet_string_.findAll(sub_sequence.alphabet_string_); }

  // Count all alphabet symbols in sequence.
  [[nodiscard]] std::vector<std::pair<typename Alphabet::Alphabet, size_t>> countSymbols() const { return alphabet_string_.countSymbols(); }

  // ReturnType the count of ordered alphabet pairs (generally used for counting CpG islands).
  [[nodiscard]] size_t countTwoSymbols(typename Alphabet::Alphabet first_symbol, typename Alphabet::Alphabet second_symbol) const { return alphabet_string_.countTwoSymbols(first_symbol, second_symbol); }

  // find the longest common prefix.
  [[nodiscard]] size_t commonPrefix(const AlphabetSequence& cmp_sequence) const { return alphabet_string_.commonPrefix(cmp_sequence.alphabet_string_); }
  // find the longest common suffix.
  [[nodiscard]] size_t commonSuffix(const AlphabetSequence& cmp_sequence) const { return alphabet_string_.commonSuffix(cmp_sequence.alphabet_string_); }
  // Extract the sequence prefix, suffix and mid-section. Generally used with commonPrefix() and commonSuffix() above.
  [[nodiscard]] AlphabetSequence removePrefix(size_t prefix_size) const { return AlphabetSequence(alphabet_string_.removePrefix(prefix_size)); }
  [[nodiscard]] AlphabetSequence removeSuffix(size_t suffix_size) const { return AlphabetSequence(alphabet_string_.removeSuffix(suffix_size)); }
  [[nodiscard]] AlphabetSequence removePrefixSuffix(size_t prefix_size, size_t suffix_size) const { return AlphabetSequence(alphabet_string_.removePrefixSuffix(prefix_size, suffix_size)); }

protected:

  AlphabetString<Alphabet> alphabet_string_;

  // Letter offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool modifyLetter(ContigOffset_t sequence_offset, typename Alphabet::Alphabet letter);
  // Delete offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool deleteOffset(ContigOffset_t delete_offset, ContigSize_t delete_size);
  // Insert offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool insertOffset(ContigOffset_t insert_offset, const AlphabetSequence& inserted_sequence);
  // Returns bool false if offset and/or size are out of bounds.
  [[nodiscard]] bool getSubsequence(ContigOffset_t substring_offset, ContigSize_t substring_size, AlphabetSequence& sub_sequence) const;
  // Equality of sequence.
  [[nodiscard]] bool equal(const AlphabetSequence& cmp_seq) const { return alphabet_string_ == cmp_seq.alphabet_string_; }
  // Equality of sub-sequence.
  [[nodiscard]] bool compareSubSequence(ContigOffset_t offset, const AlphabetSequence& sub_sequence) const {

    return alphabet_string_.compareSubString(offset, sub_sequence);

  }
  [[nodiscard]] bool compareletter(ContigOffset_t offset, typename Alphabet::Alphabet letter) const {

    return alphabet_string_.compareLetter(offset, letter);

  }

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
                                                AlphabetSequence& sub_sequence) const {

  if ((substring_offset + substring_size) > length()) {

    ExecEnv::log().error("Attempt to get substring past the end a sequence string, offset: {}, size: {}, sequence size: {}",
                         substring_offset, substring_size, length());
    return false;

  }

  if (not alphabet_string_.substring(substring_offset, substring_size, sub_sequence.alphabet_string_)) {

    ExecEnv::log().error("Problem generating a subsequence , offset: {}, subsequence size: {}, from sequence size: {}",
                         substring_offset, substring_size, alphabet_string_.length());
    return false;

  }

  return true;

}



}   // end namespace

#endif //KGL_SEQUENCE_VIRTUAL_H
