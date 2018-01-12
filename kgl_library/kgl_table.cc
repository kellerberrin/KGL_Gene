//
// Created by kellerberrin on 22/10/17.
//

#include "kgl_table.h"

namespace kgl = kellerberrin::genome;



// The NCBI Amino Acid translation tables. table should be in the interval [1, 31].
bool kgl::AminoTranslationTable::setTranslationTable(const std::string& table_name) {

  std::string uc_table_name = table_name;
  std::transform(uc_table_name.begin(), uc_table_name.end(), uc_table_name.begin(), ::toupper);
  bool table_found = false;

  for (size_t idx = 0; idx < Tables::TABLEARRAYSIZE; ++idx) {

    if (Tables::TABLEARRAY[idx]->table_name == uc_table_name) {

      amino_table_rows_ = *(Tables::TABLEARRAY[idx]);
      table_found = true;

    }

  }

  if(not table_found) {

    ExecEnv::log().warn("Amino translation table: {} not found, Standard table 1 used", table_name);
    amino_table_rows_ = *(Tables::STANDARDTABLE);

  }

  return table_found;

}


bool kgl::AminoTranslationTable::isStartAmino(AminoAcid::Alphabet amino) const {

  bool found = false;
  for (size_t index = 0; index < Tables::AMINO_TABLE_SIZE; ++index) {

    if (AminoAcid::convertToChar(amino) == amino_table_rows_.amino_table[index].amino_acid
        and amino_table_rows_.amino_table[index].start == AminoAcid::START_CODON) {

      found = true;
      break;

    }

  }

  return found;

}


size_t kgl::AminoTranslationTable::index(const Codon& codon) const {

  size_t table_index = (CodingDNA5::nucleotideToColumn(codon[0]) * Tables::CODING_NUCLEOTIDE_1) +
                       (CodingDNA5::nucleotideToColumn(codon[1]) * Tables::CODING_NUCLEOTIDE_2) +
                       CodingDNA5::nucleotideToColumn(codon[2]);

  if (table_index >= Tables::AMINO_TABLE_SIZE) {

    ExecEnv::log().error("Bad Amino Table Index: {}, codon: {}", table_index, codon.getSequenceAsString());
    table_index = amino_table_rows_.stop_codon_index;

  }

  return table_index;

}
