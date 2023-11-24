//
// Created by kellerberrin on 20/11/23.
//

#ifndef KGL_ALPHABET_SEQUENCE_H
#define KGL_ALPHABET_SEQUENCE_H


#include "kgl_alphabet_string.h"
#include "kgl_sequence_virtual.h"
#include "kel_interval_unsigned.h"

#include <string>
#include <set>


namespace kellerberrin::genome {   //  organization level namespace



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
  void clear() { alphabet_string_.clear(); }

  [[nodiscard]] ContigSize_t length() const { return alphabet_string_.length(); }
  [[nodiscard]] OpenRightUnsigned interval() const { return {0, alphabet_string_.length() }; }

  // Assumes that sequence alphabets are 1 byte (char) and map onto ascii char types. Avoids a sequence byte copy and conversion.
  [[nodiscard]] std::string_view getStringView() const override { return std::string_view{alphabet_string_.c_str(), alphabet_string_.length()}; }

  [[nodiscard]] const AlphabetString<Alphabet>& getAlphabetString() const { return alphabet_string_; }

  // Check for memory corruption.
  [[nodiscard]] bool verifySequence() const { return alphabet_string_.verifyString(); }

  // A hash value of the sequence.
  [[nodiscard]] size_t hashSequence() const { return alphabet_string_.hashString(); }

  // Insert offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool append(const AlphabetSequence& inserted_sequence);

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
  // Equality of sub-sequence.
  [[nodiscard]] bool compareSubSequence(ContigOffset_t offset, const AlphabetSequence& sub_sequence) const {

    return alphabet_string_.compareSubString(offset, sub_sequence.alphabet_string_);

  }
  [[nodiscard]] bool compareletter(ContigOffset_t offset, typename Alphabet::Alphabet letter) const {

    return alphabet_string_.compareLetter(offset, letter);

  }

  // Ptr to the base of the alphabet string. Used to initialize a SequenceView object.
  [[nodiscard]] const Alphabet::Alphabet* data() const { return alphabet_string_.data(); }

protected:

  AlphabetString<Alphabet> alphabet_string_;

  // Letter offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool modifyLetter(ContigOffset_t sequence_offset, typename Alphabet::Alphabet letter);
  // Delete offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool deleteOffset(ContigOffset_t delete_offset, ContigSize_t delete_size);
  // Insert offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool insertOffset(ContigOffset_t insert_offset, const AlphabetSequence& inserted_sequence);
  // Returns bool false if offset and/or size are out of bounds.
  [[nodiscard]] std::optional<AlphabetSequence> getSubsequence(const OpenRightUnsigned& sub_interval) const;
  // Equality of sequence.
  [[nodiscard]] bool equal(const AlphabetSequence& cmp_seq) const { return alphabet_string_ == cmp_seq.alphabet_string_; }

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
bool AlphabetSequence<Alphabet>::append(const AlphabetSequence& inserted_sequence) {

  if (not alphabet_string_.append(inserted_sequence.alphabet_string_)) {

    ExecEnv::log().error("Problem appending a sub-sequence, interval: {}, appended interval: {}",
                         interval().toString(), inserted_sequence.interval().toString());
    return false;

  }

  return true;

}


template<typename Alphabet>
bool AlphabetSequence<Alphabet>::insertOffset(ContigOffset_t insert_offset, const AlphabetSequence& inserted_sequence) {

  if (insert_offset < length()) {

    if (not alphabet_string_.insert(insert_offset, inserted_sequence.alphabet_string_)) {

      ExecEnv::log().error("Problem inserting a sub-sequence , insert offset: {}, insert sub-sequence: {}",
                           insert_offset, inserted_sequence.getStringView());
      return false;

    }

  } else if (insert_offset == length()) {

    if (not alphabet_string_.append(inserted_sequence.alphabet_string_)) {

      ExecEnv::log().error("Problem appending a sub-sequence , insert offset: {}, insert sub-sequence: {}",
                           insert_offset, inserted_sequence.getStringView());
      return false;

    }

  } else {

    ExecEnv::log().error("Attempt to insert past the end a sequence string length:{}, insert offset: {}, insert sub-sequence: {}",
                         length(), insert_offset, inserted_sequence.getStringView());
    return false;

  }


  return true;

}


template<typename Alphabet>
std::optional<AlphabetSequence<Alphabet>> AlphabetSequence<Alphabet>::getSubsequence(const OpenRightUnsigned& sub_interval) const {

  if (not interval().containsInterval(sub_interval)) {

    ExecEnv::log().warn("Sub interval: {} not contained in interval: {}.", sub_interval.toString(), interval().toString());
    return std::nullopt;

  }

  AlphabetSequence<Alphabet> sub_sequence;
  if (not alphabet_string_.substring(sub_interval.lower(), sub_interval.size(), sub_sequence.alphabet_string_)) {

    ExecEnv::log().error("Problem generating a subsequence , offset: {}, subsequence size: {}, from sequence size: {}",
                         sub_interval.lower(),sub_interval.size(), alphabet_string_.length());
    return std::nullopt;

  }

  return sub_sequence;

}


}   // end namespace



#endif //KGL_ALPHABET_SEQUENCE_H
