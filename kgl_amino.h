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

//phenylalanine - f
  aminoAcid[3][3][3] = 'f';
  aminoAcid[3][3][1] = 'F';
//Leucine - L
  aminoAcid[3][3][0] = 'L';
  aminoAcid[3][3][2] = 'L';
//Serine - S
  aminoAcid[3][1][3] = 'S';
  aminoAcid[3][1][1] = 'S';
  aminoAcid[3][1][0] = 'S';
  aminoAcid[3][1][2] = 'S';
//tyrosine - Y
  aminoAcid[3][0][3] = 'Y';
  aminoAcid[3][0][1] = 'Y';
//stop codon
  aminoAcid[3][0][0] = '-';  // OCHRE
  aminoAcid[3][0][2] = '-';  // AMBER
//cysteine - C
  aminoAcid[3][2][3] = 'C';
  aminoAcid[3][2][1] = 'C';
//stop codon
  aminoAcid[3][2][0] = '-';  // OPAL
//tryptophan - W
  aminoAcid[3][2][2] = 'W';
//leucine - L
  aminoAcid[1][3][3] = 'L';
  aminoAcid[1][3][1] = 'L';
  aminoAcid[1][3][0] = 'L';
  aminoAcid[1][3][2] = 'L';
//proline - P
  aminoAcid[1][1][3] = 'P';
  aminoAcid[1][1][1] = 'P';
  aminoAcid[1][1][0] = 'P';
  aminoAcid[1][1][2] = 'P';
//histidine - H
  aminoAcid[1][0][3] = 'H';
  aminoAcid[1][0][1] = 'H';
//glutamine - Q
  aminoAcid[1][0][0] = 'Q';
  aminoAcid[1][0][2] = 'Q';
//arginine - R
  aminoAcid[1][2][3] = 'R';
  aminoAcid[1][2][1] = 'R';
  aminoAcid[1][2][0] = 'R';
  aminoAcid[1][2][2] = 'R';
//isoleucine - I
  aminoAcid[0][3][3] = 'I';
  aminoAcid[0][3][1] = 'I';
  aminoAcid[0][3][0] = 'I';
//methionine(start codon) - M
  aminoAcid[0][3][2] = 'M';
//threonine -T
  aminoAcid[0][1][3] = 'T';
  aminoAcid[0][1][1] = 'T';
  aminoAcid[0][1][0] = 'T';
  aminoAcid[0][1][2] = 'T';
//asparagine - N
  aminoAcid[0][0][3] = 'N';
  aminoAcid[0][0][1] = 'N';
//lysine - K
  aminoAcid[0][0][0] = 'K';
  aminoAcid[0][0][2] - 'K';
//serine - S
  aminoAcid[0][2][3] = 'S';
  aminoAcid[0][2][1] = 'S';
//arginine - R
  aminoAcid[0][2][0] = 'R';
  aminoAcid[0][2][2] = 'R';
//valine - V
  aminoAcid[2][3][3] = 'V';
  aminoAcid[2][3][1] = 'V';
  aminoAcid[2][3][0] = 'V';
  aminoAcid[2][3][2] = 'V';
//alanine - A
  aminoAcid[2][1][3] = 'A';
  aminoAcid[2][1][1] = 'A';
  aminoAcid[2][1][0] = 'A';
  aminoAcid[2][1][2] = 'A';
//aspartic acid - D
  aminoAcid[2][0][3] = 'D';
  aminoAcid[2][0][1] = 'D';
//glutamic acid - E
  aminoAcid[2][0][0] = 'E';
  aminoAcid[2][0][2] = 'E';
//glycine - G
  aminoAcid[2][2][3] = 'G';
  aminoAcid[2][2][1] = 'G';
  aminoAcid[2][2][0] = 'G';
  aminoAcid[2][2][2] = 'G';

};

}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_AMINO_H
