//
// Created by kellerberrin on 17/11/17.
//

#ifndef KGL_ALPHABET_DNA5_H
#define KGL_ALPHABET_DNA5_H



#include <cstdint>
#include <string>
#include "kgl_logging.h"
#include "kgl_genome_types.h"
#include "kgl_exec_env.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Standard 5 letter nucleotide.

class DNA5 {

public:


  DNA5() = delete; // Singleton
  ~DNA5() = delete;

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

  // The Alphabet convertChar(char) function must be defined -see kgl_alphabet_string.h
  static Alphabet convertChar(char chr_base);

  // Find complementary bases.
  static Alphabet complementNucleotide(Alphabet nucleotide);

  // Convert a base to an array offset.
  static ContigOffset_t nucleotideToColumn(Alphabet nucleotide);

  // Return nucleotide as a char.
  static char convertToChar(Alphabet nucleotide) { return static_cast<char>(nucleotide); }

};




}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_ALPHABET_DNA5_H
