//
// Created by kellerberrin on 22/10/17.
//

#include "kgl_table.h"

namespace kgl = kellerberrin::genome;



// The NCBI Amino Acid translation tables. table should be in the interval [1, 31].
bool kgl::AminoTranslationTable::setTranslationTable(const std::string& table_name) {

  std::string uc_table_name = table_name;
  std::transform(uc_table_name.begin(), uc_table_name.end(), uc_table_name.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

  // Note that all table names are unique, so first match == last match.
  const auto table_iter = std::ranges::find_if(Tables::TABLEARRAY,
                                               [&uc_table_name](const TranslationTable* table) { return table->table_name == uc_table_name; });

  if (table_iter == Tables::TABLEARRAY.end()) {

    ExecEnv::log().warn("Amino translation table: {} not found, Standard table 1 used", table_name);
    amino_table_rows_ = *(Tables::STANDARDTABLE);
    return false;

  }

  amino_table_rows_ = **table_iter;
  return true;

}


bool kgl::AminoTranslationTable::isStartAmino(AminoAcid::Alphabet amino) const {

  // Both conjuncts are load-bearing: the amino_acid match alone would return true
  // for every amino that has some start codon (M and L under table 1).
  const char amino_char = AminoAcid::convertToChar(amino);
  const auto amino_table = std::span(amino_table_rows_.amino_table, Tables::AMINO_TABLE_SIZE);
  return std::ranges::any_of(amino_table, [amino_char](const AminoTableColumn& col) {

    return col.amino_acid == amino_char and col.start == AminoAcid::START_CODON;

  });

}


size_t kgl::AminoTranslationTable::index(const Codon& codon) const {

  if (codon.containsBaseN()) {

    return CONTAINS_BASE_N;

  }

  size_t table_index = (CodingDNA5::symbolToColumn(codon[0]) * Tables::CODING_NUCLEOTIDE_1) +
                       (CodingDNA5::symbolToColumn(codon[1]) * Tables::CODING_NUCLEOTIDE_2) +
                       CodingDNA5::symbolToColumn(codon[2]);

  if (table_index >= Tables::AMINO_TABLE_SIZE) {

    ExecEnv::log().error("Out of Range Amino Table Index: {}, codon: {}, (int) base1: {}, (int) base2: {}, (int) base3: {}",
                         table_index, codon.getSequenceAsString(),
                         static_cast<int>(codon[0]), static_cast<int>(codon[1]), static_cast<int>(codon[2]));
    table_index = CONTAINS_BASE_N;

  }

  return table_index;

}


kgl::AminoAcid::Alphabet kgl::AminoTranslationTable::getAmino(const Codon& codon) const {

  const size_t table_index = index(codon);

  if (table_index == CONTAINS_BASE_N) {

    return AminoAcid::AMINO_UNKNOWN;

  }

  return AminoAcid::convertChar(amino_table_rows_.amino_table[table_index].amino_acid);

}


bool kgl::AminoTranslationTable::isStopCodon(const Codon& codon) const {

  const size_t table_index = index(codon);

  if (table_index == CONTAINS_BASE_N) {

    return false;

  }

  return (amino_table_rows_.amino_table[table_index].start == AminoAcid::STOP_CODON);

}


bool kgl::AminoTranslationTable::isStartCodon(const Codon& codon) const {

  const size_t table_index = index(codon);

  if (table_index == CONTAINS_BASE_N) {

    return false;

  }

  return amino_table_rows_.amino_table[table_index].start == AminoAcid::START_CODON;

}