//
// Created by kellerberrin on 15/11/17.
//

#ifndef KGL_ALPHABET_STRING_H
#define KGL_ALPHABET_STRING_H



#include <string>
#include <algorithm>
#include "kgl_genome_types.h"
#include "kgl_alphabet_dna5.h"
#include "kgl_alphabet_amino.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Template class implements string functionality for the Nucleotide and Amino Acid alphabets.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename Alphabet>
class AlphabetString {

public:

  explicit AlphabetString() = default;
  explicit AlphabetString(const std::string& alphabet_str) { convertFromCharString(alphabet_str); }
  ~AlphabetString() = default;

  // Iterators to access the underlying std::basic_string
  using const_iterator = typename std::basic_string<typename Alphabet::Alphabet>::const_iterator;
  using const_reverse_iterator = typename std::basic_string<typename Alphabet::Alphabet>::const_reverse_iterator;
  using value_type = typename Alphabet::Alphabet;

  const_iterator begin() const { return base_string_.begin(); }
  const_reverse_iterator rbegin() const { return base_string_.rbegin(); }
  void push_back(typename Alphabet::Alphabet nucleotide) { base_string_.push_back(nucleotide); }
  void pop_back() { base_string_.pop_back(); }

  AlphabetString& operator=(const AlphabetString& copy) = default;

  ContigSize_t length() const { return base_string_.length(); }
  bool empty() const { return base_string_.empty(); }

  void reserve(ContigSize_t string_size) { base_string_.reserve(string_size); }

  bool erase(ContigOffset_t offset, ContigSize_t size) {

    try {

      base_string_.erase(offset, size);
      return true;

    }
    catch(...) {

      return false;

    }

  }

  bool insert(ContigOffset_t offset, const AlphabetString& sub_string) {

    try {

      base_string_.insert(offset, sub_string.base_string_);
      return true;

    }
    catch(...) {

      return false;

    }

  }

  typename Alphabet::Alphabet operator[] (ContigOffset_t& offset) const { return base_string_[offset]; }

  void modifyLetter(ContigOffset_t &offset, typename Alphabet::Alphabet letter) { base_string_[offset] = letter; }

  std::string str() const { return convertToCharString(); }

private:

  std::basic_string<typename Alphabet::Alphabet> base_string_;

  std::string convertToCharString() const;
  void convertFromCharString(const std::string &alphabet_str);


};


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




}
}


#endif //KGL_ALPHABET_STRING_H
