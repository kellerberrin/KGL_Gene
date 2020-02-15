//
// Created by kellerberrin on 10/01/18.
//

#ifndef KGL_TABLE_ORGANISM_H
#define KGL_TABLE_ORGANISM_H


#include "kgl_genome_types.h"
#include "kgl_table_impl.h"


namespace kellerberrin::genome {   //  organization level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines Organism Specific DNA/RNA to Amino Acid translation tables.
// Plasmodium Falciparum codon table downloaded from http://plasmodb.org (same as the standard table)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class OrganismTables {

public:

  OrganismTables() = delete; // Singleton
  ~OrganismTables() = delete;


  constexpr static const int PF_CODON_OFFSET = 48;
  constexpr static const char *PF_TABLE_NAME = "P_FALCIPARUM";
  constexpr static const char *PF_TABLE_DESC = "Organism Plasmodium Falciparum codon table downloaded from http://plasmodb.org";
  constexpr static const AminoTableColumn PFTranslationTable[Tables::AMINO_TABLE_SIZE]
  { { 'K', '-', 'A', 'A', 'A'},
  { 'N', '-', 'A', 'A', 'C'},
  { 'K', '-', 'A', 'A', 'G'},
  { 'N', '-', 'A', 'A', 'T'},
  { 'T', '-', 'A', 'C', 'A'},
  { 'T', '-', 'A', 'C', 'C'},
  { 'T', '-', 'A', 'C', 'G'},
  { 'T', '-', 'A', 'C', 'T'},
  { 'R', '-', 'A', 'G', 'A'},
  { 'S', '-', 'A', 'G', 'C'},
  { 'R', '-', 'A', 'G', 'G'},
  { 'S', '-', 'A', 'G', 'T'},
  { 'I', '-', 'A', 'T', 'A'},
  { 'I', '-', 'A', 'T', 'C'},
  { 'M', 'M', 'A', 'T', 'G'},
  { 'I', '-', 'A', 'T', 'T'},
  { 'Q', '-', 'C', 'A', 'A'},
  { 'H', '-', 'C', 'A', 'C'},
  { 'Q', '-', 'C', 'A', 'G'},
  { 'H', '-', 'C', 'A', 'T'},
  { 'P', '-', 'C', 'C', 'A'},
  { 'P', '-', 'C', 'C', 'C'},
  { 'P', '-', 'C', 'C', 'G'},
  { 'P', '-', 'C', 'C', 'T'},
  { 'R', '-', 'C', 'G', 'A'},
  { 'R', '-', 'C', 'G', 'C'},
  { 'R', '-', 'C', 'G', 'G'},
  { 'R', '-', 'C', 'G', 'T'},
  { 'L', '-', 'C', 'T', 'A'},
  { 'L', '-', 'C', 'T', 'C'},
  { 'L', '-', 'C', 'T', 'G'},
  { 'L', '-', 'C', 'T', 'T'},
  { 'E', '-', 'G', 'A', 'A'},
  { 'D', '-', 'G', 'A', 'C'},
  { 'E', '-', 'G', 'A', 'G'},
  { 'D', '-', 'G', 'A', 'T'},
  { 'A', '-', 'G', 'C', 'A'},
  { 'A', '-', 'G', 'C', 'C'},
  { 'A', '-', 'G', 'C', 'G'},
  { 'A', '-', 'G', 'C', 'T'},
  { 'G', '-', 'G', 'G', 'A'},
  { 'G', '-', 'G', 'G', 'C'},
  { 'G', '-', 'G', 'G', 'G'},
  { 'G', '-', 'G', 'G', 'T'},
  { 'V', '-', 'G', 'T', 'A'},
  { 'V', '-', 'G', 'T', 'C'},
  { 'V', '-', 'G', 'T', 'G'},
  { 'V', '-', 'G', 'T', 'T'},
  { '*', '*', 'T', 'A', 'A'},
  { 'Y', '-', 'T', 'A', 'C'},
  { '*', '*', 'T', 'A', 'G'},
  { 'Y', '-', 'T', 'A', 'T'},
  { 'S', '-', 'T', 'C', 'A'},
  { 'S', '-', 'T', 'C', 'C'},
  { 'S', '-', 'T', 'C', 'G'},
  { 'S', '-', 'T', 'C', 'T'},
  { '*', '*', 'T', 'G', 'A'},
  { 'C', '-', 'T', 'G', 'C'},
  { 'W', '-', 'T', 'G', 'G'},
  { 'C', '-', 'T', 'G', 'T'},
  { 'L', '-', 'T', 'T', 'A'},
  { 'F', '-', 'T', 'T', 'C'},
  { 'L', '-', 'T', 'T', 'G'},
  { 'F', '-', 'T', 'T', 'T'} };

  // Check the array size at compile time.
  static_assert( sizeof(PFTranslationTable)/sizeof(AminoTableColumn) == Tables::AMINO_TABLE_SIZE
  , "Error - the standard amino acid translation table should have 64 elements");

  // Define the standard translation table.
  constexpr static const TranslationTable P_FALCIPARUM{PFTranslationTable, PF_TABLE_NAME, PF_TABLE_DESC, PF_CODON_OFFSET};


};  // OrganismTable



}   // end namespace



#endif //KGL_TABLE_ORGANISM_H
