//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_ALPHABET_AMINO_H
#define KGL_ALPHABET_AMINO_H


#include "kel_exec_env.h"
#include "kgl_genome_types.h"

namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines the standard Amino Acids.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AminoAcid {

public:

  AminoAcid() = delete; // Singleton
  ~AminoAcid() = delete;

  // The 20 'natural amino acids
  inline static constexpr Amino_t PHENYLALANINE = 'F';
  inline static constexpr ContigOffset_t PHENYLALANINE_OFFSET = 0;
  inline static constexpr Amino_t LEUCINE = 'L';
  inline static constexpr ContigOffset_t LEUCINE_OFFSET = 1;
  inline static constexpr Amino_t SERINE = 'S';
  inline static constexpr ContigOffset_t SERINE_OFFSET = 2;
  inline static constexpr Amino_t TYROSINE = 'Y';
  inline static constexpr ContigOffset_t TYROSINE_OFFSET = 3;
  inline static constexpr Amino_t CYSTEINE = 'C';
  inline static constexpr ContigOffset_t CYSTEINE_OFFSET = 4;
  inline static constexpr Amino_t TRYPTOPHAN = 'W';
  inline static constexpr ContigOffset_t TRYPTOPHAN_OFFSET = 5;
  inline static constexpr Amino_t PROLINE = 'P';
  inline static constexpr ContigOffset_t PROLINE_OFFSET = 6;
  inline static constexpr Amino_t HISTIDINE = 'H';
  inline static constexpr ContigOffset_t HISTIDINE_OFFSET = 7;
  inline static constexpr Amino_t GLUTAMINE = 'Q';
  inline static constexpr ContigOffset_t GLUTAMINE_OFFSET = 8;
  inline static constexpr Amino_t ARGININE = 'R';
  inline static constexpr ContigOffset_t ARGININE_OFFSET = 9;
  inline static constexpr Amino_t ISOLEUCINE = 'I';
  inline static constexpr ContigOffset_t ISOLEUCINE_OFFSET = 10;
  inline static constexpr Amino_t METHIONINE = 'M';
  inline static constexpr ContigOffset_t METHIONINE_OFFSET = 11;
  inline static constexpr Amino_t START_CODON = 'M'; // An alias for METHIONINE
  inline static constexpr ContigOffset_t START_CODON_OFFSET = 11;
  inline static constexpr Amino_t THREONINE = 'T';
  inline static constexpr ContigOffset_t THREONINE_OFFSET = 12;
  inline static constexpr Amino_t ASPARAGINE = 'N';
  inline static constexpr ContigOffset_t ASPARAGINE_OFFSET = 13;
  inline static constexpr Amino_t LYSINE = 'K';
  inline static constexpr ContigOffset_t LYSINE_OFFSET = 14;
  inline static constexpr Amino_t VALINE = 'V';
  inline static constexpr ContigOffset_t VALINE_OFFSET = 15;
  inline static constexpr Amino_t ALANINE = 'A';
  inline static constexpr ContigOffset_t ALANINE_OFFSET = 16;
  inline static constexpr Amino_t ASPARTIC = 'D';
  inline static constexpr ContigOffset_t ASPARTIC_OFFSET = 17;
  inline static constexpr Amino_t GLUTAMIC = 'E';
  inline static constexpr ContigOffset_t GLUTAMIC_OFFSET = 18;
  inline static constexpr Amino_t GLYCINE = 'G';
  inline static constexpr ContigOffset_t GLYCINE_OFFSET = 19;
  // The additional two amino acids encoded using stop codons by some species.
  inline static constexpr Amino_t SELENOCYSTEINE = 'U';
  inline static constexpr Amino_t PYRROLYSINE = 'O';
  // The three stop codons.
  inline static constexpr Amino_t STOP_CODON = '*';
  inline static constexpr ContigOffset_t STOP_CODON_OFFSET = 20;
  inline static constexpr Amino_t STOP_AMBER_TAG = '*';
  inline static constexpr Amino_t STOP_OCHRE_TAA = '*';
  inline static constexpr Amino_t STOP_OPAL_TGA = '*';
  // Null Amino Type used in the translation tables.
  inline static constexpr Amino_t NULL_AMINO = '-';
  // The special unknown amino acid generated when the DNA5 codon
  // contains the unknown base 'N'.
  inline static constexpr Amino_t UNKNOWN_AMINO = 'Z';
  inline static constexpr ContigOffset_t UNKNOWN_AMINO_OFFSET = 21;

  // The Alphabet enum type must be defined -see kgl_alphabet_string.h
  enum class Alphabet : Amino_t
  { F = PHENYLALANINE,
    L = LEUCINE,
    S = SERINE,
    Y = TYROSINE,
    C = CYSTEINE,
    W = TRYPTOPHAN,
    P = PROLINE,
    H = HISTIDINE,
    Q = GLUTAMINE,
    R = ARGININE,
    I = ISOLEUCINE,
    M = METHIONINE,
    T = THREONINE,
    N = ASPARAGINE,
    K = LYSINE,
    V = VALINE,
    A = ALANINE,
    D = ASPARTIC,
    E = GLUTAMIC,
    G = GLYCINE,
    // Rare - The additional two amino acids encoded as stop codons by some species.
    U = SELENOCYSTEINE,
    O = PYRROLYSINE,
    // Represents all stop codons, including the two above.
    _ = STOP_CODON,   // Note '*' is an illegal symbol in an enum.
    // The special unknown amino acid generated when the CodingDNA5 codon
    // contains the unknown base 'N'.
    Z = UNKNOWN_AMINO};

  inline static constexpr Alphabet AMINO_STOP = Alphabet::_;
  inline static constexpr Alphabet AMINO_UNKNOWN = Alphabet::Z;

  // Return a vector of all valid alphabet values (only the 20 'natural'
  // and the 'unknown'  amino acids are returned).
  [[nodiscard]] static const std::vector<Alphabet>& enumerateAlphabet();

  // Checks for possible memory corruption.
  [[nodiscard]] static bool validAlphabet(Alphabet amino);

  // Convert an amino into an array offset.
  // Should be in the same ascending order as the vector returned by enumerateAlphabet();
  // Only the 20 'natural' AAs and the unknown AA are indexed.
  [[nodiscard]] static ContigOffset_t symbolToColumn(Alphabet amino);

  // The Alphabet convertChar(char) function must be defined -see kgl_alphabet_string.h
  [[nodiscard]] static Alphabet convertChar(char chr_aa);

  // Return amino acid as a char.
  [[nodiscard]] static char convertToChar(Alphabet amino) { return static_cast<char>(amino); }

};


}   // end namespace


#endif //KGL_ALPHABET_AMINO_H
