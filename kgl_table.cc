//
// Created by kellerberrin on 22/10/17.
//

#include "kgl_table.h"

namespace kgl = kellerberrin::genome;


// The NCBI Amino Acid translation tables. table should be in the interval [1, 31].
bool kgl::AminoTranslationTable::setTranslationTable(size_t table) {

  bool result = false;

  if(table > Tables::TABLE_ARRAY_SIZE or table == 0) {

    ExecEnv::log().warn("Amino translation table: {} not implemented, Standard table 1 used", table);
    table = Tables::STANDARD_TABLE_1;
    result = false;

  }
  if (Tables::TABLEARRAY[table-1] == nullptr) {

    ExecEnv::log().warn("Amino translation table: {} not implemented, Standard table 1 used", table);
    table = Tables::STANDARD_TABLE_1;
    result = false;

  }

  amino_table_rows_ = *(Tables::TABLEARRAY[table-1]);

  return result;

}


size_t kgl::AminoTranslationTable::index(const AminoAcidTypes::Codon& Codon) {

  size_t table_index = (NucleotideColumn_DNA5::nucleotideToColumn(Codon.bases[0]) * Tables::CODING_NUCLEOTIDE_1) +
                       (NucleotideColumn_DNA5::nucleotideToColumn(Codon.bases[1]) * Tables::CODING_NUCLEOTIDE_2) +
                       NucleotideColumn_DNA5::nucleotideToColumn(Codon.bases[2]);

  if (table_index >= Tables::AMINO_TABLE_SIZE) {

    ExecEnv::log().error("Bad Amino Table Index: {}, base1: {} base2: {}, base3: {}",
                         table_index, Codon.bases[0], Codon.bases[1], Codon.bases[2]);
    table_index = amino_table_rows_.stop_codon_index;

  }

  return table_index;

}
