//
// Created by kellerberrin on 18/10/17.
//

#ifndef KGL_AMINO_H
#define KGL_AMINO_H

#include "kgl_exec_env.h"
#include "kgl_nucleotide.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

class AminoAcidTypes {

public:

  using AminoType = Amino_22_t;

  AminoAcidTypes() = delete; // Singleton
  ~AminoAcidTypes() = delete;

  static constexpr size_t CODON_SIZE = 3;
  static constexpr size_t CODING_NUCLEOTIDES = 4;

  struct Codon {

    const NucleotideColumn_DNA5::NucleotideType* bases;

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

};



// Amino acid Translation tables
// Found at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

class StandardAminoTranslationTable {

public:

  StandardAminoTranslationTable() = delete; // Singleton
  ~StandardAminoTranslationTable() = delete;

  static inline AminoAcidTypes::AminoType getAmino(const AminoAcidTypes::Codon& Codon) {

    return amino_table_rows_[index(Codon)].amino_acid;

  }

  static inline bool isStopCodon(const AminoAcidTypes::Codon& Codon) {

    if (Codon.bases[0] == 'T' and Codon.bases[1] == 'A' and Codon.bases[2] == 'G') {

      ExecEnv::log().info("Check for stop codon with 'TAG', table index: {}, Amino: {}, bases: {}{}{}, flag: {}",
                          index(Codon),
                          amino_table_rows_[index(Codon)].amino_acid,
                          amino_table_rows_[index(Codon)].base1,
                          amino_table_rows_[index(Codon)].base2,
                          amino_table_rows_[index(Codon)].base3,
                          (amino_table_rows_[index(Codon)].start == AminoAcidTypes::STOP_CODON));

    }

    return (amino_table_rows_[index(Codon)].start == AminoAcidTypes::STOP_CODON);

  }

  static inline bool isStartCodon(const AminoAcidTypes::Codon& Codon) {

    return amino_table_rows_[index(Codon)].start == AminoAcidTypes::START_CODON;

  }


private:

  static inline size_t index(const AminoAcidTypes::Codon& Codon) {


    size_t table_index = (NucleotideColumn_DNA5::nucleotideToColumn(Codon.bases[0]) * CODING_NUCLEOTIDE_1) +
    (NucleotideColumn_DNA5::nucleotideToColumn(Codon.bases[1]) * CODING_NUCLEOTIDE_2) +
    NucleotideColumn_DNA5::nucleotideToColumn(Codon.bases[2]);

    if (table_index >= AMINO_TABLE_SIZE) {

      ExecEnv::log().error("Bad Amino Table Index: {}, base1: {} base2: {}, base3: {}",
                           table_index, Codon.bases[0], Codon.bases[1], Codon.bases[2]);
      table_index = STOP_CODON_OFFSET;
    }

    return table_index;

  }

  constexpr static int AMINO_TABLE_SIZE = 64;
  constexpr static int CODING_NUCLEOTIDE_2 = 4;
  constexpr static int CODING_NUCLEOTIDE_1 = 16;
  constexpr static int STOP_CODON_OFFSET = 48;

  struct AminoTableColumn {
    AminoAcidTypes::AminoType amino_acid;
    AminoAcidTypes::AminoType start;
    NucleotideColumn_DNA5::NucleotideType base1;
    NucleotideColumn_DNA5::NucleotideType base2;
    NucleotideColumn_DNA5::NucleotideType base3;
  };

  constexpr static AminoTableColumn amino_table_rows_[AMINO_TABLE_SIZE] =
  {{ 'K', '-', 'A', 'A', 'A' } ,
   { 'N', '-', 'A', 'A', 'C' } ,
   { 'K', '-', 'A', 'A', 'G' } ,
   { 'N', '-', 'A', 'A', 'T' } ,
   { 'T', '-', 'A', 'C', 'A' } ,
   { 'T', '-', 'A', 'C', 'C' } ,
   { 'T', '-', 'A', 'C', 'G' } ,
   { 'T', '-', 'A', 'C', 'T' } ,
   { 'R', '-', 'A', 'G', 'A' } ,
   { 'S', '-', 'A', 'G', 'C' } ,
   { 'R', '-', 'A', 'G', 'G' } ,
   { 'S', '-', 'A', 'G', 'T' } ,
   { 'I', '-', 'A', 'T', 'A' } ,
   { 'I', '-', 'A', 'T', 'C' } ,
   { 'M', 'M', 'A', 'T', 'G' } ,
   { 'I', '-', 'A', 'T', 'T' } ,
   { 'Q', '-', 'C', 'A', 'A' } ,
   { 'H', '-', 'C', 'A', 'C' } ,
   { 'Q', '-', 'C', 'A', 'G' } ,
   { 'H', '-', 'C', 'A', 'T' } ,
   { 'P', '-', 'C', 'C', 'A' } ,
   { 'P', '-', 'C', 'C', 'C' } ,
   { 'P', '-', 'C', 'C', 'G' } ,
   { 'P', '-', 'C', 'C', 'T' } ,
   { 'R', '-', 'C', 'G', 'A' } ,
   { 'R', '-', 'C', 'G', 'C' } ,
   { 'R', '-', 'C', 'G', 'G' } ,
   { 'R', '-', 'C', 'G', 'T' } ,
   { 'L', '-', 'C', 'T', 'A' } ,
   { 'L', '-', 'C', 'T', 'C' } ,
   { 'L', 'M', 'C', 'T', 'G' } ,
   { 'L', '-', 'C', 'T', 'T' } ,
   { 'E', '-', 'G', 'A', 'A' } ,
   { 'D', '-', 'G', 'A', 'C' } ,
   { 'E', '-', 'G', 'A', 'G' } ,
   { 'D', '-', 'G', 'A', 'T' } ,
   { 'A', '-', 'G', 'C', 'A' } ,
   { 'A', '-', 'G', 'C', 'C' } ,
   { 'A', '-', 'G', 'C', 'G' } ,
   { 'A', '-', 'G', 'C', 'T' } ,
   { 'G', '-', 'G', 'G', 'A' } ,
   { 'G', '-', 'G', 'G', 'C' } ,
   { 'G', '-', 'G', 'G', 'G' } ,
   { 'G', '-', 'G', 'G', 'T' } ,
   { 'V', '-', 'G', 'T', 'A' } ,
   { 'V', '-', 'G', 'T', 'C' } ,
   { 'V', '-', 'G', 'T', 'G' } ,
   { 'V', '-', 'G', 'T', 'T' } ,
   { '*', '*', 'T', 'A', 'A' } ,
   { 'Y', '-', 'T', 'A', 'C' } ,
   { '*', '*', 'T', 'A', 'G' } ,
   { 'Y', '-', 'T', 'A', 'T' } ,
   { 'S', '-', 'T', 'C', 'A' } ,
   { 'S', '-', 'T', 'C', 'C' } ,
   { 'S', '-', 'T', 'C', 'G' } ,
   { 'S', '-', 'T', 'C', 'T' } ,
   { '*', '*', 'T', 'G', 'A' } ,
   { 'C', '-', 'T', 'G', 'C' } ,
   { 'W', '-', 'T', 'G', 'G' } ,
   { 'C', '-', 'T', 'G', 'T' } ,
   { 'L', '-', 'T', 'T', 'A' } ,
   { 'F', '-', 'T', 'T', 'C' } ,
   { 'L', 'M', 'T', 'T', 'G' } ,
   { 'F', '-', 'T', 'T', 'T' }};

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_AMINO_H
