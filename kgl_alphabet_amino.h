//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_ALPHABET_AMINO_H
#define KGL_ALPHABET_AMINO_H


#include "kgl_exec_env.h"
#include "kgl_genome_types.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines the standard Amino Acids.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AminoAcid {

public:

  AminoAcid() = delete; // Singleton
  ~AminoAcid() = delete;

  // The 20 'natural amino acids
  static constexpr Amino_t PHENYLALANINE = 'F';
  static constexpr Amino_t LEUCINE = 'L';
  static constexpr Amino_t SERINE = 'S';
  static constexpr Amino_t TYROSINE = 'Y';
  static constexpr Amino_t CYSTEINE = 'C';
  static constexpr Amino_t TRYPTOPHAN = 'W';
  static constexpr Amino_t PROLINE = 'P';
  static constexpr Amino_t HISTIDINE = 'H';
  static constexpr Amino_t GLUTAMINE = 'Q';
  static constexpr Amino_t ARGININE = 'R';
  static constexpr Amino_t ISOLEUCINE = 'I';
  static constexpr Amino_t METHIONINE = 'M';
  static constexpr Amino_t START_CODON = 'M';
  static constexpr Amino_t THREONINE = 'T';
  static constexpr Amino_t ASPARAGINE = 'N';
  static constexpr Amino_t LYSINE = 'K';
  static constexpr Amino_t VALINE = 'V';
  static constexpr Amino_t ALANINE = 'A';
  static constexpr Amino_t ASPARTIC = 'D';
  static constexpr Amino_t GLUTAMIC = 'E';
  static constexpr Amino_t GLYCINE = 'G';
  // The additional two amino acids encoded using stop codons by some species.
  static constexpr Amino_t SELENOCYSTEINE = 'U';
  static constexpr Amino_t PYRROLYSINE = 'O';
  // The three stop codons.
  static constexpr Amino_t STOP_CODON = '*';
  static constexpr Amino_t STOP_AMBER_TAG = '*';
  static constexpr Amino_t STOP_OCHRE_TAA = '*';
  static constexpr Amino_t STOP_OPAL_TGA = '*';
  // Null Amino Type used in the translation tables.
  static constexpr Amino_t NULL_AMINO = '-';
  // The special unknown amino acid generated when the DNA5 codon
  // contains the unknown base 'N'.
  static constexpr Amino_t UNKNOWN_AMINO = '?';

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
    // Rare - The additional two amino acids encoded using stop codons by some species.
    U = SELENOCYSTEINE,
    O = PYRROLYSINE,
    // The three stop codons.
    _ = STOP_CODON,
    // The special unknown amino acid generated when the DNA5 codon
    // contains the unknown base 'N'.
    Z = UNKNOWN_AMINO};

  static constexpr Alphabet AMINO_STOP = Alphabet::_;
  static constexpr Alphabet AMINO_UNKNOWN = Alphabet::Z;

  // The Alphabet convertChar(char) function must be defined -see kgl_alphabet_string.h
  static Alphabet convertChar(char chr_base);

  // Return amin acid as a char.
  static char convertToChar(Alphabet nucleotide) { return static_cast<char>(nucleotide); }

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_ALPHABET_AMINO_H
