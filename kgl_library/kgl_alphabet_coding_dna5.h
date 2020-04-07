//
// Created by kellerberrin on 25/11/17.
//

#ifndef KGL_ALPHABET_CODING_DNA5_H
#define KGL_ALPHABET_CODING_DNA5_H




#include <cstdint>
#include <string>
#include "kel_logging.h"
#include "kgl_genome_types.h"
#include "kel_exec_env.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a semantic class that defines strand converted DNA sequences.
// Sequences that contain strings of this class have been STRANDED and can be used to generate amino acids.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CodingDNA5 {

public:


  CodingDNA5() = delete; // Singleton
  ~CodingDNA5() = delete;

  inline static constexpr ContigOffset_t NUCLEOTIDE_COLUMNS = 5;

  inline static constexpr Nucleotide_t A_NUCLEOTIDE = 'A';
  inline static constexpr ContigOffset_t A_NUCLEOTIDE_OFFSET = 0;
  inline static constexpr Nucleotide_t C_NUCLEOTIDE = 'C';
  inline static constexpr ContigOffset_t C_NUCLEOTIDE_OFFSET = 1;
  inline static constexpr Nucleotide_t G_NUCLEOTIDE = 'G';
  inline static constexpr ContigOffset_t G_NUCLEOTIDE_OFFSET = 2;
  inline static constexpr Nucleotide_t U_NUCLEOTIDE = 'U';
  inline static constexpr ContigOffset_t U_NUCLEOTIDE_OFFSET = 3;
  inline static constexpr Nucleotide_t T_NUCLEOTIDE = 'T';
  inline static constexpr ContigOffset_t T_NUCLEOTIDE_OFFSET = 3;
  inline static constexpr Nucleotide_t N_NUCLEOTIDE = 'N';
  inline static constexpr ContigOffset_t N_NUCLEOTIDE_OFFSET = 4;

  // The Alphabet enum type must be defined -see kgl_alphabet_string.h
  enum class Alphabet : Nucleotide_t {
    A = A_NUCLEOTIDE,
    C = C_NUCLEOTIDE,
    G = G_NUCLEOTIDE,
    T = T_NUCLEOTIDE,
    N = N_NUCLEOTIDE
  };

  // Return a vector of all valid alphabet values.
  [[nodiscard]] static std::vector<Alphabet> enumerateAlphabet();

  // Checks for possible memory corruption.
  [[nodiscard]] static bool validAlphabet(Alphabet nucleotide);

  // Convert a base to an array offset.
  [[nodiscard]] static ContigOffset_t nucleotideToColumn(Alphabet nucleotide);

  // Return nucleotide as a char.
  [[nodiscard]] static char convertToChar(Alphabet nucleotide) { return static_cast<char>(nucleotide); }



};


}   // end namespace



#endif //KGL_ALPHABET_CODING_DNA5_H
