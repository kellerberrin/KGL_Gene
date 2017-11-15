//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_ALPHABET_AMINO_H
#define KGL_ALPHABET_AMINO_H


#include "kgl_exec_env.h"
#include "kgl_alphabet_base.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines the standard Amino Acids. Used as a template class with the Amino Acid sequence objects.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AminoAcid {

public:

  AminoAcid() = delete; // Singleton
  ~AminoAcid() = delete;

  static constexpr size_t CODING_NUCLEOTIDES = 4;

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
  static constexpr Amino_t LYSINE = 'L';
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


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_ALPHABET_AMINO_H
