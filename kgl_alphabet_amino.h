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

class AminoAcidTypes {

public:

  using AminoType = Amino_22_t;

  AminoAcidTypes() = delete; // Singleton
  ~AminoAcidTypes() = delete;

  static constexpr size_t CODON_SIZE = 3;
  static constexpr size_t CODING_NUCLEOTIDES = 4;

  struct Codon {

    const Nucleotide_ExtendedDNA5* bases;

  };

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
  static constexpr AminoType START_CODON = 'M';
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
  static constexpr AminoType STOP_CODON = '*';
  static constexpr AminoType STOP_AMBER_TAG = '*';
  static constexpr AminoType STOP_OCHRE_TAA = '*';
  static constexpr AminoType STOP_OPAL_TGA = '*';
  // Null Amino Type used in the translation tables.
  static constexpr AminoType NULL_AMINO = '-';
  // Lower case amino acids
  static constexpr AminoType PHENYLALANINE_LC = 'f';
  static constexpr AminoType LEUCINE_LC = 'l';
  static constexpr AminoType SERINE_LC = 's';
  static constexpr AminoType TYROSINE_LC = 'y';
  static constexpr AminoType CYSTEINE_LC = 'c';
  static constexpr AminoType TRYPTOPHAN_LC = 'w';
  static constexpr AminoType PROLINE_LC = 'p';
  static constexpr AminoType HISTIDINE_LC = 'h';
  static constexpr AminoType GLUTAMINE_LC = 'q';
  static constexpr AminoType ARGININE_LC = 'r';
  static constexpr AminoType ISOLEUCINE_LC = 'i';
  static constexpr AminoType METHIONINE_LC = 'm';
  static constexpr AminoType THREONINE_LC = 't';
  static constexpr AminoType ASPARAGINE_LC = 'n';
  static constexpr AminoType LYSINE_LC = 'l';
  static constexpr AminoType VALINE_LC = 'v';
  static constexpr AminoType ALANINE_LC = 'a';
  static constexpr AminoType ASPARTIC_LC = 'd';
  static constexpr AminoType GLUTAMIC_LC = 'e';
  static constexpr AminoType GLYCINE_LC = 'g';
  // The additional two amino acids encoded using stop codons by some species.
  static constexpr AminoType SELENOCYSTEINE_LC = 'u';
  static constexpr AminoType PYRROLYSINE_LC = 'o';


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_ALPHABET_AMINO_H
