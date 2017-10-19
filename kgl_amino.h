//
// Created by kellerberrin on 18/10/17.
//

#ifndef KGL_AMINO_H
#define KGL_AMINO_H

#include "kgl_nucleotide.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

class AminoColumn_64 {

public:

  using AminoType = Amino_22_t;

  AminoColumn_64() = delete; // Singleton
  ~AminoColumn_64() = delete;

//RNA codon to amino acid mapping
//A = 0, C = 1, G = 2, U = 3

  // The 20 'natural amino acids
  static constexpr AminoType PHENYLALANINE = 'F';
  static constexpr AminoType LEUCINE = 'L';
  static constexpr AminoType SERINE = 'S';
  static constexpr AminoType TYROSINE = 'Y';
  static constexpr AminoType CYSTEINE = 'C';
  static constexpr AminoType TRYPTOPHAN = 'W';
  static constexpr AminoType PROLINE = 'P';
  static constexpr AminoType HISTIDINE = 'H';
  static constexpr AminoType GLUTAMINE = 'Q';
  static constexpr AminoType ARGININE = 'R';
  static constexpr AminoType ISOLEUCINE = 'I';
  static constexpr AminoType METHIONINE = 'M';
  static constexpr AminoType THREONINE = 'T';
  static constexpr AminoType ASPARAGINE = 'N';
  static constexpr AminoType LYSINE = 'L';
  static constexpr AminoType VALINE = 'V';
  static constexpr AminoType ALANINE = 'A';
  static constexpr AminoType ASPARTIC = 'D';
  static constexpr AminoType GLUTAMIC = 'E';
  static constexpr AminoType GLYCINE = 'G';
  // The additional two amino acids encoded using stop codons by some species.
  static constexpr AminoType SELENOCYSTEINE = 'U';
  static constexpr AminoType PYRROLYSINE = 'O';
  // The three stop codons.
  static constexpr AminoType STOP_AMBER_TAG = '*';
  static constexpr AminoType STOP_OCHRE_TAA = '*';
  static constexpr AminoType STOP_OPAL_TGA = '*';

  Nucleotide_DNA5_t aminoAcid[4][4][4];

};

}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_AMINO_H
