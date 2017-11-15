//
// Created by kellerberrin on 15/11/17.
//

#ifndef KGL_ALPHABET_STRING_H
#define KGL_ALPHABET_STRING_H



#include <string>
#include <algorithm>
#include "kgl_genome_types.h"
#include "kgl_alphabet_base.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Template class implements string functionality for the Nucleotide and Amino Acid alphabets.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename Alphabet>
class StringAlphabet {

public:

  explicit StringAlphabet(const std::string& alphabet_str) { convertFromString(alphabet_str); }
  explicit StringAlphabet(std::basic_string<typename Alphabet::Alphabet> base_string) : base_string_(std::move(base_string)) {}
  ~StringAlphabet() = default;

  StringAlphabet& operator=(const StringAlphabet& copy) = default;

  StringAlphabet substr(ContigOffset_t offset, ContigSize_t size)
  {

    return StringAlphabet(base_string_.substr(offset, size));

  }

  ContigSize_t length() const { return base_string_.length(); }

  std::string str() const { return convertToString(); }

private:

  std::basic_string<typename Alphabet::Alphabet> base_string_;

  std::string convertToString() const;
  void convertFromString(const std::string& alphabet_str);


};


template<typename Alphabet>
std::string StringAlphabet<Alphabet>::convertToString() const {

  std::string base_string;
  base_string.reserve(length());

  auto convert_base = [](typename Alphabet::Alphabet base) { return static_cast<char>(base); };
  std::transform(base_string_.begin(), base_string_.end(), std::back_inserter(base_string), convert_base);

  return base_string;

}


template<typename Alphabet>
void StringAlphabet<Alphabet>::convertFromString(const std::string &alphabet_str) {

  base_string_.reserve(alphabet_str.length());
  std::transform(alphabet_str.begin(), alphabet_str.end(), std::back_inserter(base_string_), Alphabet::convertChar);

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The alphabet strings are defined here.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The standard 5 nucleotide DNA/RNA A, C, G, T/U, N
using StringDNA5 = StringAlphabet<DNA5>;

// The extended DNA alphabet to include deletions and insertions.
// Note that X = Delete, E = A insert, F = C insert, I = G insert, J = T insert and K = N insert.
using StringExtendDNA5 = StringAlphabet<ExtendDNA5>;




}
}


#endif //KGL_ALPHABET_STRING_H
