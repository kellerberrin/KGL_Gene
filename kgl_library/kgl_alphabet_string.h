//
// Created by kellerberrin on 15/11/17.
//

#ifndef KGL_ALPHABET_STRING_H
#define KGL_ALPHABET_STRING_H



#include <string>
#include <algorithm>
#include <functional>
#include "kgl_genome_types.h"
#include "kgl_alphabet_dna5.h"
#include "kgl_alphabet_amino.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Template class implements string functionality for the Nucleotide and Amino Acid alphabets.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename Alphabet>
class AlphabetString {

public:

  explicit AlphabetString() = default;
  AlphabetString(AlphabetString<Alphabet>&& alphabet_string) noexcept : base_string_(std::move(alphabet_string.base_string_)) {}
  AlphabetString(const AlphabetString<Alphabet>& alphabet_string) : base_string_(alphabet_string.base_string_) {}
  explicit AlphabetString(const std::string& alphabet_str) { convertFromCharString(alphabet_str); }
  ~AlphabetString() = default;

  // Iterators to access the underlying std::basic_string
  using const_iterator = typename std::basic_string<typename Alphabet::Alphabet>::const_iterator;
  using const_reverse_iterator = typename std::basic_string<typename Alphabet::Alphabet>::const_reverse_iterator;
  using value_type = typename Alphabet::Alphabet;

  [[nodiscard]] const_iterator begin() const { return base_string_.begin(); }
  [[nodiscard]] const_reverse_iterator rbegin() const { return base_string_.rbegin(); }
  [[nodiscard]] const_iterator end() const { return base_string_.end(); }
  [[nodiscard]] const_reverse_iterator rend() const { return base_string_.rend(); }
  void push_back(typename Alphabet::Alphabet nucleotide) { base_string_.push_back(nucleotide); }
  void pop_back() { base_string_.pop_back(); }

  AlphabetString& operator=(const AlphabetString& copy) = default;

  [[nodiscard]] ContigSize_t length() const { return base_string_.length(); }
  [[nodiscard]] bool empty() const { return base_string_.empty(); }

  void reserve(ContigSize_t string_size) { base_string_.reserve(string_size); }

  [[nodiscard]] bool erase(ContigOffset_t offset, ContigSize_t size) {

    try {

      base_string_.erase(offset, size);
      return true;

    }
    catch(...) {

      return false;

    }

  }

  void assignString(const std::string &alphabet_str) { convertFromCharString(alphabet_str); }

  [[nodiscard]] bool insert(ContigOffset_t offset, const AlphabetString& sub_string) {

    try {

      base_string_.insert(offset, sub_string.base_string_);
      return true;

    }
    catch(...) {

      return false;

    }

  }

  [[nodiscard]] bool append(const AlphabetString& sub_string) {

    try {

      base_string_.append(sub_string.base_string_);
      return true;

    }
    catch(...) {

      return false;

    }

  }


  [[nodiscard]] bool substring(ContigOffset_t offset, ContigSize_t size, AlphabetString& substring) const {

    try {

      substring.base_string_ = base_string_.substr(offset, size);
      return true;

    }
    catch(...) {

      return false;

    }

  }


  typename Alphabet::Alphabet operator[] (ContigOffset_t& offset) const { return base_string_[offset]; }

  void modifyLetter(ContigOffset_t &offset, typename Alphabet::Alphabet letter) { base_string_[offset] = letter; }

  [[nodiscard]] std::string str() const { return convertToCharString(); }

  bool operator==(const AlphabetString& compare_string) const { return (base_string_ == compare_string.base_string_); }

  std::vector<ContigOffset_t> findAll(const AlphabetString& sub_string) const {

    std::vector<ContigOffset_t> offset_vector;
    size_t offset = base_string_.find(sub_string.base_string_);

    while (offset != std::basic_string<Alphabet>::npos) {

      offset_vector.push_back(static_cast<ContigOffset_t>(offset));
      offset = base_string_.find(sub_string.base_string_, offset + 1);

    }

    return offset_vector;

  }

  [[nodiscard]] bool verifyString() const;

  [[nodiscard]] size_t hashString() const;

private:

  std::basic_string<typename Alphabet::Alphabet> base_string_;

  [[nodiscard]] std::string convertToCharString() const;
  void convertFromCharString(const std::string &alphabet_str);


};


template<typename Alphabet>
size_t AlphabetString<Alphabet>::hashString() const {

  return std::hash<std::string>{}(str());

}



template<typename Alphabet>
bool AlphabetString<Alphabet>::verifyString() const {

  for (ContigSize_t idx = 0; idx < length(); ++idx) {

    bool compare = Alphabet::validAlphabet(base_string_.at(idx));

    if (not compare) {

      ExecEnv::log().error("AlphabetString::verifyString(), invalid Alphabet value (int): {} found at index: {}", static_cast<size_t>(base_string_.at(idx)), idx);
      return false;

    }

  }

  return true;

}


template<typename Alphabet>
std::string AlphabetString<Alphabet>::convertToCharString() const {

  std::string base_string;
  base_string.reserve(length());

  auto convert_base = [](typename Alphabet::Alphabet base) { return static_cast<char>(base); };
  std::transform(base_string_.begin(), base_string_.end(), std::back_inserter(base_string), convert_base);

  return base_string;

}


template<typename Alphabet>
void AlphabetString<Alphabet>::convertFromCharString(const std::string &alphabet_str) {

  base_string_.reserve(alphabet_str.length());
  std::transform(alphabet_str.begin(), alphabet_str.end(), std::back_inserter(base_string_), Alphabet::convertChar);

}




}  // end namespace


#endif //KGL_ALPHABET_STRING_H
