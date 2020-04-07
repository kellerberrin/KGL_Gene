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
  inline static constexpr Amino_t LEUCINE = 'L';
  inline static constexpr Amino_t SERINE = 'S';
  inline static constexpr Amino_t TYROSINE = 'Y';
  inline static constexpr Amino_t CYSTEINE = 'C';
  inline static constexpr Amino_t TRYPTOPHAN = 'W';
  inline static constexpr Amino_t PROLINE = 'P';
  inline static constexpr Amino_t HISTIDINE = 'H';
  inline static constexpr Amino_t GLUTAMINE = 'Q';
  inline static constexpr Amino_t ARGININE = 'R';
  inline static constexpr Amino_t ISOLEUCINE = 'I';
  inline static constexpr Amino_t METHIONINE = 'M';
  inline static constexpr Amino_t START_CODON = 'M';
  inline static constexpr Amino_t THREONINE = 'T';
  inline static constexpr Amino_t ASPARAGINE = 'N';
  inline static constexpr Amino_t LYSINE = 'K';
  inline static constexpr Amino_t VALINE = 'V';
  inline static constexpr Amino_t ALANINE = 'A';
  inline static constexpr Amino_t ASPARTIC = 'D';
  inline static constexpr Amino_t GLUTAMIC = 'E';
  inline static constexpr Amino_t GLYCINE = 'G';
  // The additional two amino acids encoded using stop codons by some species.
  inline static constexpr Amino_t SELENOCYSTEINE = 'U';
  inline static constexpr Amino_t PYRROLYSINE = 'O';
  // The three stop codons.
  inline static constexpr Amino_t STOP_CODON = '*';
  inline static constexpr Amino_t STOP_AMBER_TAG = '*';
  inline static constexpr Amino_t STOP_OCHRE_TAA = '*';
  inline static constexpr Amino_t STOP_OPAL_TGA = '*';
  // Null Amino Type used in the translation tables.
  inline static constexpr Amino_t NULL_AMINO = '-';
  // The special unknown amino acid generated when the DNA5 codon
  // contains the unknown base 'N'.
  inline static constexpr Amino_t UNKNOWN_AMINO = '?';

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
    _ = STOP_CODON,
    // The special unknown amino acid generated when the CodingDNA5 codon
    // contains the unknown base 'N'.
    Z = UNKNOWN_AMINO};

  inline static constexpr Alphabet AMINO_STOP = Alphabet::_;
  inline static constexpr Alphabet AMINO_UNKNOWN = Alphabet::Z;

  // Return a vector of all valid alphabet values (only the 20 'natural'
  // and the 'unknown'  amino acids are returned).
  [[nodiscard]] static std::vector<Alphabet> enumerateAlphabet();

  // Checks for possible memory corruption.
  [[nodiscard]] static bool validAlphabet(Alphabet nucleotide);

  // The Alphabet convertChar(char) function must be defined -see kgl_alphabet_string.h
  [[nodiscard]] static Alphabet convertChar(char chr_base);

  // Return amino acid as a char.
  [[nodiscard]] static char convertToChar(Alphabet amino) { return static_cast<char>(amino); }

};


}   // end namespace


#endif //KGL_ALPHABET_AMINO_H
