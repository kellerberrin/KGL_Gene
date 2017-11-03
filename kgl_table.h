//
// Created by kellerberrin on 22/10/17.
//

#ifndef KGL_TABLE_H
#define KGL_TABLE_H


#include "kgl_table_ncbi.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coding Table Class
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Amino acid Translation tables
// Found at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
class AminoTranslationTable {

public:

  explicit AminoTranslationTable() { setTranslationTable(Tables::STANDARD_TABLE_1); }
  ~AminoTranslationTable() = default;

  std::string TableName() { return amino_table_rows_.table_name; }

  // The NCBI Amino Acid translation tables. table should be in the interval [1, 31].
  bool setTranslationTable(size_t table);

  AminoAcidTypes::AminoType getAmino(const AminoAcidTypes::Codon& Codon) {

    return amino_table_rows_.amino_table[index(Codon)].amino_acid;

  }

  AminoAcidTypes::AminoType getStopAmino() {

    return amino_table_rows_.amino_table[amino_table_rows_.stop_codon_index].amino_acid;

  }

  bool isStopCodon(const AminoAcidTypes::Codon& Codon) {

    return (amino_table_rows_.amino_table[index(Codon)].start == AminoAcidTypes::STOP_CODON);

  }

  bool isStartCodon(const AminoAcidTypes::Codon& Codon) {

    return amino_table_rows_.amino_table[index(Codon)].start == AminoAcidTypes::START_CODON;

  }

private:

  TranslationTable amino_table_rows_;

  size_t index(const AminoAcidTypes::Codon& Codon);

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_TABLE_H
