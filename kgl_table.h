//
// Created by kellerberrin on 22/10/17.
//

#ifndef KGL_TABLE_H
#define KGL_TABLE_H


#include "kgl_nucleotide.h"
#include "kgl_amino.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



struct AminoTableColumn {

  AminoAcidTypes::AminoType amino_acid;
  AminoAcidTypes::AminoType start;
  NucleotideColumn_DNA5::NucleotideType base1;
  NucleotideColumn_DNA5::NucleotideType base2;
  NucleotideColumn_DNA5::NucleotideType base3;

};

struct TranslationTable {

  const AminoTableColumn* amino_table;
  const char* table_name;
  const size_t stop_codon_index;

};


struct Tables {

  constexpr static const int AMINO_TABLE_SIZE = 64;

  constexpr static const int TABLE_1_STOP = 48;
  constexpr static const char *TABLE_1_NAME = "The Standard Code (table 1)";
  constexpr static const AminoTableColumn TABLE_1_TABLE[AMINO_TABLE_SIZE]
  {{'K', '-', 'A', 'A', 'A'},
   {'N', '-', 'A', 'A', 'C'},
   {'K', '-', 'A', 'A', 'G'},
   {'N', '-', 'A', 'A', 'T'},
   {'T', '-', 'A', 'C', 'A'},
   {'T', '-', 'A', 'C', 'C'},
   {'T', '-', 'A', 'C', 'G'},
   {'T', '-', 'A', 'C', 'T'},
   {'R', '-', 'A', 'G', 'A'},
   {'S', '-', 'A', 'G', 'C'},
   {'R', '-', 'A', 'G', 'G'},
   {'S', '-', 'A', 'G', 'T'},
   {'I', '-', 'A', 'T', 'A'},
   {'I', '-', 'A', 'T', 'C'},
   {'M', 'M', 'A', 'T', 'G'},
   {'I', '-', 'A', 'T', 'T'},
   {'Q', '-', 'C', 'A', 'A'},
   {'H', '-', 'C', 'A', 'C'},
   {'Q', '-', 'C', 'A', 'G'},
   {'H', '-', 'C', 'A', 'T'},
   {'P', '-', 'C', 'C', 'A'},
   {'P', '-', 'C', 'C', 'C'},
   {'P', '-', 'C', 'C', 'G'},
   {'P', '-', 'C', 'C', 'T'},
   {'R', '-', 'C', 'G', 'A'},
   {'R', '-', 'C', 'G', 'C'},
   {'R', '-', 'C', 'G', 'G'},
   {'R', '-', 'C', 'G', 'T'},
   {'L', '-', 'C', 'T', 'A'},
   {'L', '-', 'C', 'T', 'C'},
   {'L', 'M', 'C', 'T', 'G'},
   {'L', '-', 'C', 'T', 'T'},
   {'E', '-', 'G', 'A', 'A'},
   {'D', '-', 'G', 'A', 'C'},
   {'E', '-', 'G', 'A', 'G'},
   {'D', '-', 'G', 'A', 'T'},
   {'A', '-', 'G', 'C', 'A'},
   {'A', '-', 'G', 'C', 'C'},
   {'A', '-', 'G', 'C', 'G'},
   {'A', '-', 'G', 'C', 'T'},
   {'G', '-', 'G', 'G', 'A'},
   {'G', '-', 'G', 'G', 'C'},
   {'G', '-', 'G', 'G', 'G'},
   {'G', '-', 'G', 'G', 'T'},
   {'V', '-', 'G', 'T', 'A'},
   {'V', '-', 'G', 'T', 'C'},
   {'V', '-', 'G', 'T', 'G'},
   {'V', '-', 'G', 'T', 'T'},
   {'*', '*', 'T', 'A', 'A'},
   {'Y', '-', 'T', 'A', 'C'},
   {'*', '*', 'T', 'A', 'G'},
   {'Y', '-', 'T', 'A', 'T'},
   {'S', '-', 'T', 'C', 'A'},
   {'S', '-', 'T', 'C', 'C'},
   {'S', '-', 'T', 'C', 'G'},
   {'S', '-', 'T', 'C', 'T'},
   {'*', '*', 'T', 'G', 'A'},
   {'C', '-', 'T', 'G', 'C'},
   {'W', '-', 'T', 'G', 'G'},
   {'C', '-', 'T', 'G', 'T'},
   {'L', '-', 'T', 'T', 'A'},
   {'F', '-', 'T', 'T', 'C'},
   {'L', 'M', 'T', 'T', 'G'},
   {'F', '-', 'T', 'T', 'T'}};

  constexpr static const TranslationTable TABLE_1{TABLE_1_TABLE, TABLE_1_NAME, TABLE_1_STOP};

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_TABLE_H
