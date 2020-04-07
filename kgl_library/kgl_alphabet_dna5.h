//
// Created by kellerberrin on 17/11/17.
//

#ifndef KGL_ALPHABET_DNA5_H
#define KGL_ALPHABET_DNA5_H



#include <cstdint>
#include <string>
#include "kel_logging.h"
#include "kgl_genome_types.h"
#include "kel_exec_env.h"
#include "kgl_alphabet_coding_dna5.h"


namespace kellerberrin::genome {   //  organization::project level namespace


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

  // The Alphabet convertChar(char) function must be defined -see kgl_alphabet_string.h
  [[nodiscard]] static Alphabet convertChar(char chr_base);

  // Find complementary bases.
  [[nodiscard]] static CodingDNA5::Alphabet complementNucleotide(Alphabet nucleotide);


  // Convert to CodingDNA5 without complementary base conversion.
  // Warning - assumes that CodingDNA5::Alphabet and DNA5::Alphabet nucleotides have the same enum values.
  [[nodiscard]] static CodingDNA5::Alphabet convertToCodingDNA5(Alphabet nucleotide) { return static_cast<CodingDNA5::Alphabet>(nucleotide); }

  // Convert from coding DNA5.
  // Warning - assumes that CodingDNA5::Alphabet and DNA5::Alphabet nucleotides have the same enum values.
  [[nodiscard]] static DNA5::Alphabet convertFromCodingDNA5(CodingDNA5::Alphabet nucleotide) { return static_cast<DNA5::Alphabet>(nucleotide); }

  // Find complementary bases and convert to DNA5.
  [[nodiscard]] static DNA5::Alphabet convertComplementNucleotide(CodingDNA5::Alphabet nucleotide) {

    return convertFromCodingDNA5(complementNucleotide(convertFromCodingDNA5(nucleotide)));

  }

  // Convert a base to an array offset.
  [[nodiscard]] static ContigOffset_t nucleotideToColumn(Alphabet nucleotide);

  // Return nucleotide as a char.
  [[nodiscard]] static char convertToChar(Alphabet nucleotide) { return static_cast<char>(nucleotide); }

  // Converts an array offset into a base.
  [[nodiscard]] static Alphabet offsetToNucleotide(ContigOffset_t offset) {

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

  }

};




}   // end namespace



#endif //KGL_ALPHABET_DNA5_H
