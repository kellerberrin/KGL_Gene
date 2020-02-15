//
// Created by kellerberrin on 25/11/17.
//

#ifndef KGL_ALPHABET_CODING_DNA5_H
#define KGL_ALPHABET_CODING_DNA5_H




#include <cstdint>
#include <string>
#include "kgl_logging.h"
#include "kgl_genome_types.h"
#include "kgl_exec_env.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a semantic class that defines strand converted DNA sequences.
// Sequences that contain strings of this class have been STRANDED and can be used to generate amino acids.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CodingDNA5 {

public:


  CodingDNA5() = delete; // Singleton
  ~CodingDNA5() = delete;

  static constexpr ContigOffset_t NUCLEOTIDE_COLUMNS = 5;

  static constexpr Nucleotide_t A_NUCLEOTIDE = 'A';
  static constexpr ContigOffset_t A_NUCLEOTIDE_OFFSET = 0;
  static constexpr Nucleotide_t C_NUCLEOTIDE = 'C';
  static constexpr ContigOffset_t C_NUCLEOTIDE_OFFSET = 1;
  static constexpr Nucleotide_t G_NUCLEOTIDE = 'G';
  static constexpr ContigOffset_t G_NUCLEOTIDE_OFFSET = 2;
  static constexpr Nucleotide_t U_NUCLEOTIDE = 'U';
  static constexpr ContigOffset_t U_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_t T_NUCLEOTIDE = 'T';
  static constexpr ContigOffset_t T_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_t N_NUCLEOTIDE = 'N';
  static constexpr ContigOffset_t N_NUCLEOTIDE_OFFSET = 4;


  // The Alphabet enum type must be defined -see kgl_alphabet_string.h
  enum class Alphabet : Nucleotide_t {
    A = A_NUCLEOTIDE,
    C = C_NUCLEOTIDE,
    G = G_NUCLEOTIDE,
    T = T_NUCLEOTIDE,
    N = N_NUCLEOTIDE
  };

  // Return a vector of all valid alphabet values.
  static std::vector<Alphabet> enumerateAlphabet();

  // Checks for possible memory corruption.
  static bool validAlphabet(Alphabet nucleotide);

  // Convert a base to an array offset.
  static ContigOffset_t nucleotideToColumn(Alphabet nucleotide);

  // Return nucleotide as a char.
  static char convertToChar(Alphabet nucleotide) { return static_cast<char>(nucleotide); }



};


}   // end namespace



#endif //KGL_ALPHABET_CODING_DNA5_H
