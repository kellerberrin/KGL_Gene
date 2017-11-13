//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_TABLE_NCBI_H
#define KGL_TABLE_NCBI_H


#include "kgl_alphabet_base.h"
#include "kgl_alphabet_amino.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines DNA/RNA to Amino Acid translation tables.
// These are found at the NCBI website: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
// At the current time only 5 tables are implemented (there are 31) as it is a somewhat tedious task.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct AminoTableColumn {

  AminoAcidTypes::AminoType amino_acid;
  AminoAcidTypes::AminoType start;
  Nucleotide_DNA5_t base1;
  Nucleotide_DNA5_t base2;
  Nucleotide_DNA5_t base3;

};


struct TranslationTable {

  const AminoTableColumn* amino_table;
  const char* table_name;
  size_t stop_codon_index;

};



class Tables {

public:

  Tables() = delete; // Singleton
  ~Tables() = delete;

  constexpr static const int CODING_NUCLEOTIDE_2 = 4;
  constexpr static const int CODING_NUCLEOTIDE_1 = 16;
  constexpr static const int AMINO_TABLE_SIZE = 64;
  constexpr static const int STANDARD_TABLE_1 = 1;

private:

  constexpr static const int STOP_CODON_OFFSET = 48;
  constexpr static const char *AMINO_TABLE_NAME = "The Standard Code (table 1)";
  constexpr static const AminoTableColumn StandardTranslationTable[AMINO_TABLE_SIZE]
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

  // Check the array size at compile time.
  static_assert( sizeof(StandardTranslationTable)/sizeof(AminoTableColumn) == AMINO_TABLE_SIZE
  , "Error - the standard amino acid translation table should have 64 elements");

  // Define the standard translation table.
  constexpr static const TranslationTable TABLE_1{StandardTranslationTable, AMINO_TABLE_NAME, STOP_CODON_OFFSET};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Vertebrate Mitochondrial Code (transl_table=2)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  constexpr static const int STOP_CODON_OFFSET_2 = 48;
  constexpr static const char *AMINO_TABLE_2_NAME = "The Vertebrate Mitochondrial Code (transl_table=2)";
  constexpr static const AminoTableColumn TranslationTable_2[AMINO_TABLE_SIZE]
  {{ 'K' ,'-' ,'A' ,'A' ,'A' },
   { 'N' ,'-' ,'A' ,'A' ,'C' },
   { 'K' ,'-' ,'A' ,'A' ,'G' },
   { 'N' ,'-' ,'A' ,'A' ,'T' },
   { 'T' ,'-' ,'A' ,'C' ,'A' },
   { 'T' ,'-' ,'A' ,'C' ,'C' },
   { 'T' ,'-' ,'A' ,'C' ,'G' },
   { 'T' ,'-' ,'A' ,'C' ,'T' },
   { '*' ,'*' ,'A' ,'G' ,'A' },
   { 'S' ,'-' ,'A' ,'G' ,'C' },
   { '*' ,'*' ,'A' ,'G' ,'G' },
   { 'S' ,'-' ,'A' ,'G' ,'T' },
   { 'M' ,'M' ,'A' ,'T' ,'A' },
   { 'I' ,'M' ,'A' ,'T' ,'C' },
   { 'M' ,'M' ,'A' ,'T' ,'G' },
   { 'I' ,'M' ,'A' ,'T' ,'T' },
   { 'Q' ,'-' ,'C' ,'A' ,'A' },
   { 'H' ,'-' ,'C' ,'A' ,'C' },
   { 'Q' ,'-' ,'C' ,'A' ,'G' },
   { 'H' ,'-' ,'C' ,'A' ,'T' },
   { 'P' ,'-' ,'C' ,'C' ,'A' },
   { 'P' ,'-' ,'C' ,'C' ,'C' },
   { 'P' ,'-' ,'C' ,'C' ,'G' },
   { 'P' ,'-' ,'C' ,'C' ,'T' },
   { 'R' ,'-' ,'C' ,'G' ,'A' },
   { 'R' ,'-' ,'C' ,'G' ,'C' },
   { 'R' ,'-' ,'C' ,'G' ,'G' },
   { 'R' ,'-' ,'C' ,'G' ,'T' },
   { 'L' ,'-' ,'C' ,'T' ,'A' },
   { 'L' ,'-' ,'C' ,'T' ,'C' },
   { 'L' ,'-' ,'C' ,'T' ,'G' },
   { 'L' ,'-' ,'C' ,'T' ,'T' },
   { 'E' ,'-' ,'G' ,'A' ,'A' },
   { 'D' ,'-' ,'G' ,'A' ,'C' },
   { 'E' ,'-' ,'G' ,'A' ,'G' },
   { 'D' ,'-' ,'G' ,'A' ,'T' },
   { 'A' ,'-' ,'G' ,'C' ,'A' },
   { 'A' ,'-' ,'G' ,'C' ,'C' },
   { 'A' ,'-' ,'G' ,'C' ,'G' },
   { 'A' ,'-' ,'G' ,'C' ,'T' },
   { 'G' ,'-' ,'G' ,'G' ,'A' },
   { 'G' ,'-' ,'G' ,'G' ,'C' },
   { 'G' ,'-' ,'G' ,'G' ,'G' },
   { 'G' ,'-' ,'G' ,'G' ,'T' },
   { 'V' ,'-' ,'G' ,'T' ,'A' },
   { 'V' ,'-' ,'G' ,'T' ,'C' },
   { 'V' ,'M' ,'G' ,'T' ,'G' },
   { 'V' ,'-' ,'G' ,'T' ,'T' },
   { '*' ,'*' ,'T' ,'A' ,'A' },
   { 'Y' ,'-' ,'T' ,'A' ,'C' },
   { '*' ,'*' ,'T' ,'A' ,'G' },
   { 'Y' ,'-' ,'T' ,'A' ,'T' },
   { 'S' ,'-' ,'T' ,'C' ,'A' },
   { 'S' ,'-' ,'T' ,'C' ,'C' },
   { 'S' ,'-' ,'T' ,'C' ,'G' },
   { 'S' ,'-' ,'T' ,'C' ,'T' },
   { 'W' ,'-' ,'T' ,'G' ,'A' },
   { 'C' ,'-' ,'T' ,'G' ,'C' },
   { 'W' ,'-' ,'T' ,'G' ,'G' },
   { 'C' ,'-' ,'T' ,'G' ,'T' },
   { 'L' ,'-' ,'T' ,'T' ,'A' },
   { 'F' ,'-' ,'T' ,'T' ,'C' },
   { 'L' ,'-' ,'T' ,'T' ,'G' },
   { 'F' ,'-' ,'T' ,'T' ,'T' }};

  // Check the array size at compile time.
  static_assert( sizeof(TranslationTable_2)/sizeof(AminoTableColumn) == AMINO_TABLE_SIZE
  , "Error - the standard amino acid translation table should have 64 elements");

  // Define the standard translation table.
  constexpr static const TranslationTable TABLE_2{TranslationTable_2, AMINO_TABLE_2_NAME, STOP_CODON_OFFSET_2};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Yeast Mitochondrial Code (transl_table=3)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  constexpr static const int STOP_CODON_OFFSET_3 = 48;
  constexpr static const char *AMINO_TABLE_3_NAME = "The Yeast Mitochondrial Code (transl_table=3)";
  constexpr static const AminoTableColumn TranslationTable_3[AMINO_TABLE_SIZE]
  {{ 'K' ,'-' ,'A' ,'A' ,'A' },
   { 'N' ,'-' ,'A' ,'A' ,'C' },
   { 'K' ,'-' ,'A' ,'A' ,'G' },
   { 'N' ,'-' ,'A' ,'A' ,'T' },
   { 'T' ,'-' ,'A' ,'C' ,'A' },
   { 'T' ,'-' ,'A' ,'C' ,'C' },
   { 'T' ,'-' ,'A' ,'C' ,'G' },
   { 'T' ,'-' ,'A' ,'C' ,'T' },
   { 'R' ,'-' ,'A' ,'G' ,'A' },
   { 'S' ,'-' ,'A' ,'G' ,'C' },
   { 'R' ,'-' ,'A' ,'G' ,'G' },
   { 'S' ,'-' ,'A' ,'G' ,'T' },
   { 'M' ,'M' ,'A' ,'T' ,'A' },
   { 'I' ,'-' ,'A' ,'T' ,'C' },
   { 'M' ,'M' ,'A' ,'T' ,'G' },
   { 'I' ,'-' ,'A' ,'T' ,'T' },
   { 'Q' ,'-' ,'C' ,'A' ,'A' },
   { 'H' ,'-' ,'C' ,'A' ,'C' },
   { 'Q' ,'-' ,'C' ,'A' ,'G' },
   { 'H' ,'-' ,'C' ,'A' ,'T' },
   { 'P' ,'-' ,'C' ,'C' ,'A' },
   { 'P' ,'-' ,'C' ,'C' ,'C' },
   { 'P' ,'-' ,'C' ,'C' ,'G' },
   { 'P' ,'-' ,'C' ,'C' ,'T' },
   { 'R' ,'-' ,'C' ,'G' ,'A' },
   { 'R' ,'-' ,'C' ,'G' ,'C' },
   { 'R' ,'-' ,'C' ,'G' ,'G' },
   { 'R' ,'-' ,'C' ,'G' ,'T' },
   { 'T' ,'-' ,'C' ,'T' ,'A' },
   { 'T' ,'-' ,'C' ,'T' ,'C' },
   { 'T' ,'-' ,'C' ,'T' ,'G' },
   { 'T' ,'-' ,'C' ,'T' ,'T' },
   { 'E' ,'-' ,'G' ,'A' ,'A' },
   { 'D' ,'-' ,'G' ,'A' ,'C' },
   { 'E' ,'-' ,'G' ,'A' ,'G' },
   { 'D' ,'-' ,'G' ,'A' ,'T' },
   { 'A' ,'-' ,'G' ,'C' ,'A' },
   { 'A' ,'-' ,'G' ,'C' ,'C' },
   { 'A' ,'-' ,'G' ,'C' ,'G' },
   { 'A' ,'-' ,'G' ,'C' ,'T' },
   { 'G' ,'-' ,'G' ,'G' ,'A' },
   { 'G' ,'-' ,'G' ,'G' ,'C' },
   { 'G' ,'-' ,'G' ,'G' ,'G' },
   { 'G' ,'-' ,'G' ,'G' ,'T' },
   { 'V' ,'-' ,'G' ,'T' ,'A' },
   { 'V' ,'-' ,'G' ,'T' ,'C' },
   { 'V' ,'-' ,'G' ,'T' ,'G' },
   { 'V' ,'-' ,'G' ,'T' ,'T' },
   { '*' ,'*' ,'T' ,'A' ,'A' },
   { 'Y' ,'-' ,'T' ,'A' ,'C' },
   { '*' ,'*' ,'T' ,'A' ,'G' },
   { 'Y' ,'-' ,'T' ,'A' ,'T' },
   { 'S' ,'-' ,'T' ,'C' ,'A' },
   { 'S' ,'-' ,'T' ,'C' ,'C' },
   { 'S' ,'-' ,'T' ,'C' ,'G' },
   { 'S' ,'-' ,'T' ,'C' ,'T' },
   { 'W' ,'-' ,'T' ,'G' ,'A' },
   { 'C' ,'-' ,'T' ,'G' ,'C' },
   { 'W' ,'-' ,'T' ,'G' ,'G' },
   { 'C' ,'-' ,'T' ,'G' ,'T' },
   { 'L' ,'-' ,'T' ,'T' ,'A' },
   { 'F' ,'-' ,'T' ,'T' ,'C' },
   { 'L' ,'-' ,'T' ,'T' ,'G' },
   { 'F' ,'-' ,'T' ,'T' ,'T' }};

  // Check the array size at compile time.
  static_assert( sizeof(TranslationTable_3)/sizeof(AminoTableColumn) == AMINO_TABLE_SIZE
  , "Error - the standard amino acid translation table should have 64 elements");

  // Define the standard translation table.
  constexpr static const TranslationTable TABLE_3{TranslationTable_3, AMINO_TABLE_3_NAME, STOP_CODON_OFFSET_3};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  constexpr static const int STOP_CODON_OFFSET_4 = 48;
  constexpr static const char *AMINO_TABLE_4_NAME =
  "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)";
  constexpr static const AminoTableColumn TranslationTable_4[AMINO_TABLE_SIZE]
  {{ 'K' ,'-' ,'A' ,'A' ,'A' },
   { 'N' ,'-' ,'A' ,'A' ,'C' },
   { 'K' ,'-' ,'A' ,'A' ,'G' },
   { 'N' ,'-' ,'A' ,'A' ,'T' },
   { 'T' ,'-' ,'A' ,'C' ,'A' },
   { 'T' ,'-' ,'A' ,'C' ,'C' },
   { 'T' ,'-' ,'A' ,'C' ,'G' },
   { 'T' ,'-' ,'A' ,'C' ,'T' },
   { 'R' ,'-' ,'A' ,'G' ,'A' },
   { 'S' ,'-' ,'A' ,'G' ,'C' },
   { 'R' ,'-' ,'A' ,'G' ,'G' },
   { 'S' ,'-' ,'A' ,'G' ,'T' },
   { 'I' ,'M' ,'A' ,'T' ,'A' },
   { 'I' ,'M' ,'A' ,'T' ,'C' },
   { 'M' ,'M' ,'A' ,'T' ,'G' },
   { 'I' ,'M' ,'A' ,'T' ,'T' },
   { 'Q' ,'-' ,'C' ,'A' ,'A' },
   { 'H' ,'-' ,'C' ,'A' ,'C' },
   { 'Q' ,'-' ,'C' ,'A' ,'G' },
   { 'H' ,'-' ,'C' ,'A' ,'T' },
   { 'P' ,'-' ,'C' ,'C' ,'A' },
   { 'P' ,'-' ,'C' ,'C' ,'C' },
   { 'P' ,'-' ,'C' ,'C' ,'G' },
   { 'P' ,'-' ,'C' ,'C' ,'T' },
   { 'R' ,'-' ,'C' ,'G' ,'A' },
   { 'R' ,'-' ,'C' ,'G' ,'C' },
   { 'R' ,'-' ,'C' ,'G' ,'G' },
   { 'R' ,'-' ,'C' ,'G' ,'T' },
   { 'L' ,'-' ,'C' ,'T' ,'A' },
   { 'L' ,'-' ,'C' ,'T' ,'C' },
   { 'L' ,'M' ,'C' ,'T' ,'G' },
   { 'L' ,'-' ,'C' ,'T' ,'T' },
   { 'E' ,'-' ,'G' ,'A' ,'A' },
   { 'D' ,'-' ,'G' ,'A' ,'C' },
   { 'E' ,'-' ,'G' ,'A' ,'G' },
   { 'D' ,'-' ,'G' ,'A' ,'T' },
   { 'A' ,'-' ,'G' ,'C' ,'A' },
   { 'A' ,'-' ,'G' ,'C' ,'C' },
   { 'A' ,'-' ,'G' ,'C' ,'G' },
   { 'A' ,'-' ,'G' ,'C' ,'T' },
   { 'G' ,'-' ,'G' ,'G' ,'A' },
   { 'G' ,'-' ,'G' ,'G' ,'C' },
   { 'G' ,'-' ,'G' ,'G' ,'G' },
   { 'G' ,'-' ,'G' ,'G' ,'T' },
   { 'V' ,'-' ,'G' ,'T' ,'A' },
   { 'V' ,'-' ,'G' ,'T' ,'C' },
   { 'V' ,'M' ,'G' ,'T' ,'G' },
   { 'V' ,'-' ,'G' ,'T' ,'T' },
   { '*' ,'*' ,'T' ,'A' ,'A' },
   { 'Y' ,'-' ,'T' ,'A' ,'C' },
   { '*' ,'*' ,'T' ,'A' ,'G' },
   { 'Y' ,'-' ,'T' ,'A' ,'T' },
   { 'S' ,'-' ,'T' ,'C' ,'A' },
   { 'S' ,'-' ,'T' ,'C' ,'C' },
   { 'S' ,'-' ,'T' ,'C' ,'G' },
   { 'S' ,'-' ,'T' ,'C' ,'T' },
   { 'W' ,'-' ,'T' ,'G' ,'A' },
   { 'C' ,'-' ,'T' ,'G' ,'C' },
   { 'W' ,'-' ,'T' ,'G' ,'G' },
   { 'C' ,'-' ,'T' ,'G' ,'T' },
   { 'L' ,'M' ,'T' ,'T' ,'A' },
   { 'F' ,'-' ,'T' ,'T' ,'C' },
   { 'L' ,'M' ,'T' ,'T' ,'G' },
   { 'F' ,'-' ,'T' ,'T' ,'T' }};

  // Check the array size at compile time.
  static_assert( sizeof(TranslationTable_4)/sizeof(AminoTableColumn) == AMINO_TABLE_SIZE
  , "Error - the standard amino acid translation table should have 64 elements");

  // Define the standard translation table.
  constexpr static const TranslationTable TABLE_4{TranslationTable_4, AMINO_TABLE_4_NAME, STOP_CODON_OFFSET_4};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Invertebrate Mitochondrial Code (transl_table=5)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  constexpr static const int STOP_CODON_OFFSET_5 = 48;
  constexpr static const char *AMINO_TABLE_5_NAME = "The Invertebrate Mitochondrial Code (transl_table=5)";
  constexpr static const AminoTableColumn TranslationTable_5[AMINO_TABLE_SIZE]
  {{ 'K' ,'-' ,'A' ,'A' ,'A' },
   { 'N' ,'-' ,'A' ,'A' ,'C' },
   { 'K' ,'-' ,'A' ,'A' ,'G' },
   { 'N' ,'-' ,'A' ,'A' ,'T' },
   { 'T' ,'-' ,'A' ,'C' ,'A' },
   { 'T' ,'-' ,'A' ,'C' ,'C' },
   { 'T' ,'-' ,'A' ,'C' ,'G' },
   { 'T' ,'-' ,'A' ,'C' ,'T' },
   { 'S' ,'-' ,'A' ,'G' ,'A' },
   { 'S' ,'-' ,'A' ,'G' ,'C' },
   { 'S' ,'-' ,'A' ,'G' ,'G' },
   { 'S' ,'-' ,'A' ,'G' ,'T' },
   { 'M' ,'M' ,'A' ,'T' ,'A' },
   { 'I' ,'M' ,'A' ,'T' ,'C' },
   { 'M' ,'M' ,'A' ,'T' ,'G' },
   { 'I' ,'M' ,'A' ,'T' ,'T' },
   { 'Q' ,'-' ,'C' ,'A' ,'A' },
   { 'H' ,'-' ,'C' ,'A' ,'C' },
   { 'Q' ,'-' ,'C' ,'A' ,'G' },
   { 'H' ,'-' ,'C' ,'A' ,'T' },
   { 'P' ,'-' ,'C' ,'C' ,'A' },
   { 'P' ,'-' ,'C' ,'C' ,'C' },
   { 'P' ,'-' ,'C' ,'C' ,'G' },
   { 'P' ,'-' ,'C' ,'C' ,'T' },
   { 'R' ,'-' ,'C' ,'G' ,'A' },
   { 'R' ,'-' ,'C' ,'G' ,'C' },
   { 'R' ,'-' ,'C' ,'G' ,'G' },
   { 'R' ,'-' ,'C' ,'G' ,'T' },
   { 'L' ,'-' ,'C' ,'T' ,'A' },
   { 'L' ,'-' ,'C' ,'T' ,'C' },
   { 'L' ,'-' ,'C' ,'T' ,'G' },
   { 'L' ,'-' ,'C' ,'T' ,'T' },
   { 'E' ,'-' ,'G' ,'A' ,'A' },
   { 'D' ,'-' ,'G' ,'A' ,'C' },
   { 'E' ,'-' ,'G' ,'A' ,'G' },
   { 'D' ,'-' ,'G' ,'A' ,'T' },
   { 'A' ,'-' ,'G' ,'C' ,'A' },
   { 'A' ,'-' ,'G' ,'C' ,'C' },
   { 'A' ,'-' ,'G' ,'C' ,'G' },
   { 'A' ,'-' ,'G' ,'C' ,'T' },
   { 'G' ,'-' ,'G' ,'G' ,'A' },
   { 'G' ,'-' ,'G' ,'G' ,'C' },
   { 'G' ,'-' ,'G' ,'G' ,'G' },
   { 'G' ,'-' ,'G' ,'G' ,'T' },
   { 'V' ,'-' ,'G' ,'T' ,'A' },
   { 'V' ,'-' ,'G' ,'T' ,'C' },
   { 'V' ,'M' ,'G' ,'T' ,'G' },
   { 'V' ,'-' ,'G' ,'T' ,'T' },
   { '*' ,'*' ,'T' ,'A' ,'A' },
   { 'Y' ,'-' ,'T' ,'A' ,'C' },
   { '*' ,'*' ,'T' ,'A' ,'G' },
   { 'Y' ,'-' ,'T' ,'A' ,'T' },
   { 'S' ,'-' ,'T' ,'C' ,'A' },
   { 'S' ,'-' ,'T' ,'C' ,'C' },
   { 'S' ,'-' ,'T' ,'C' ,'G' },
   { 'S' ,'-' ,'T' ,'C' ,'T' },
   { 'W' ,'-' ,'T' ,'G' ,'A' },
   { 'C' ,'-' ,'T' ,'G' ,'C' },
   { 'W' ,'-' ,'T' ,'G' ,'G' },
   { 'C' ,'-' ,'T' ,'G' ,'T' },
   { 'L' ,'-' ,'T' ,'T' ,'A' },
   { 'F' ,'-' ,'T' ,'T' ,'C' },
   { 'L' ,'M' ,'T' ,'T' ,'G' },
   { 'F' ,'-' ,'T' ,'T' ,'T' }};

  // Check the array size at compile time.
  static_assert( sizeof(TranslationTable_5)/sizeof(AminoTableColumn) == AMINO_TABLE_SIZE
  , "Error - the standard amino acid translation table should have 64 elements");

  // Define the standard translation table.
  constexpr static const TranslationTable TABLE_5{TranslationTable_5, AMINO_TABLE_5_NAME, STOP_CODON_OFFSET_5};

public:

  // An array of Amino translation tables, nullptr are ignored and default to TABLE_1 (the standard table)
  // These tables are defined on the NCBI website: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
  // There are 31 tables defined there and as time permits these can be added to the table array and will
  // be automatically available for DNA/RNA translation and Gene verification. Note that the code maintains
  // separate tables for different contigs, so for example, a different table could used for (e.g. table 2)
  // for a vertebrate mitochondrial contig.
  constexpr static const size_t TABLE_ARRAY_SIZE = 5;
  constexpr static const TranslationTable* TABLEARRAY[]{ &TABLE_1, &TABLE_2, &TABLE_3, &TABLE_4, &TABLE_5 };


};  // Table



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_TABLE_NCBI_H
