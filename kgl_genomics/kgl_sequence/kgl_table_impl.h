//
// Created by kellerberrin on 10/01/18.
//

#ifndef KGL_TABLE_IMPL_H
#define KGL_TABLE_IMPL_H


#include <vector>
#include "kgl_genome_types.h"


namespace kellerberrin::genome {   //  organization level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines DNA/RNA to Amino Acid translation tables.
// These are found at the NCBI website: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
// At the current time only 5 tables are implemented (there are 31) as it is a somewhat tedious task.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct AminoTableColumn {

  Amino_t amino_acid;
  Amino_t start;
  Nucleotide_t base1;
  Nucleotide_t base2;
  Nucleotide_t base3;

};


struct TranslationTable {

  const AminoTableColumn* amino_table;
  const char* table_name;
  const char* table_description;
  size_t stop_codon_index;

};


class Tables {

public:

  Tables(); // initialize the table vector.
  ~Tables() = default;

  constexpr static const int CODING_NUCLEOTIDE_2 = 4;
  constexpr static const int CODING_NUCLEOTIDE_1 = 16;
  constexpr static const int AMINO_TABLE_SIZE = 64;

  // These tables are defined on the NCBI website: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
  // There are 31 tables defined there and as time permits these can be added to the table array and will
  // be automatically available for DNA/RNA translation and Gene verification. Note that the kgl code maintains
  // separate tables for different contigs, so for example, a different table could used for (e.g. table 2)
  // for a vertebrate mitochondrial contig.
  static const TranslationTable* STANDARDTABLE;
  static const TranslationTable* TABLEARRAY[];
  static const size_t TABLEARRAYSIZE;

};  // Table



}   // end namespace



#endif //KGL_TABLE_IMPL_H
