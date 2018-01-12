//
// Created by kellerberrin on 22/10/17.
//

#ifndef KGL_TABLE_H
#define KGL_TABLE_H


#include "kgl_table_impl.h"
#include "kgl_sequence_codon.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coding Table Class
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Amino acid Translation tables
// Found at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
class AminoTranslationTable {

public:

  explicit AminoTranslationTable() : amino_table_rows_(*Tables::STANDARDTABLE) {}
  ~AminoTranslationTable() = default;

  std::string TableName() const { return amino_table_rows_.table_name; }

  std::string TableDescription() const { return amino_table_rows_.table_description; }

  // Get table by name.
  bool setTranslationTable(const std::string& table_name);

  AminoAcid::Alphabet getAmino(const Codon& Codon) {

    return AminoAcid::convertChar(amino_table_rows_.amino_table[index(Codon)].amino_acid);

  }

  bool isStopAmino(AminoAcid::Alphabet amino) const { return amino == AminoAcid::AMINO_STOP; }

  bool isStartAmino(AminoAcid::Alphabet amino) const;

  bool isStopCodon(const Codon& Codon) const {

    return (amino_table_rows_.amino_table[index(Codon)].start == AminoAcid::STOP_CODON);

  }

  bool isStartCodon(const Codon& Codon) const {

    return amino_table_rows_.amino_table[index(Codon)].start == AminoAcid::START_CODON;

  }

private:

  TranslationTable amino_table_rows_;

  size_t index(const Codon& Codon) const;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_TABLE_H
