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
#include "kgl_alphabet_coding_dna5.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This class defines DNA sequences.
// The sequences that contain this class have NOT been strand converted.
// Do not use this alphabet to generate amino sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
  static CodingDNA5::Alphabet complementNucleotide(Alphabet nucleotide);

  // Convert to CodingDNA5 without complementary base conversion.
  static CodingDNA5::Alphabet convertToCodingDNA5(Alphabet nucleotide) { return static_cast<CodingDNA5::Alphabet>(nucleotide); }

  // Convert from coding DNA5.
  static DNA5::Alphabet convertFromCodingDNA5(CodingDNA5::Alphabet nucleotide) { return static_cast<DNA5::Alphabet>(nucleotide); }

  // Convert a base to an array offset.
  static ContigOffset_t nucleotideToColumn(Alphabet nucleotide);

  // Return nucleotide as a char.
  static char convertToChar(Alphabet nucleotide) { return static_cast<char>(nucleotide); }

  // return bool true if G or C
  static bool isNucleotideGC(Alphabet nucleotide) { return nucleotide == Alphabet::G or nucleotide == Alphabet::C; }

  // Converts an array offset into a base.
  static Alphabet offsetToNucleotide(ContigOffset_t offset) {

    // Translate the nucleotide to an array column
    switch (offset) {

      case A_NUCLEOTIDE_OFFSET: return Alphabet::A;
      case C_NUCLEOTIDE_OFFSET: return Alphabet::C;
      case G_NUCLEOTIDE_OFFSET: return Alphabet::G;
      case T_NUCLEOTIDE_OFFSET: return Alphabet::T;
      case N_NUCLEOTIDE_OFFSET: return Alphabet::N;

      default:
        ExecEnv::log().error("DNA5::offsetToNucleotide(), Invalid Nucleotide Offset", offset);
        return Alphabet::N;
    }

    return Alphabet::N; // Never reached, to keep the compiler happy.

  }

};




}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_ALPHABET_DNA5_H
