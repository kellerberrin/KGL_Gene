//
// Created by kellerberrin on 22/10/17.
//

#ifndef KGL_TABLE_H
#define KGL_TABLE_H


#include "kgl_table_impl.h"
#include "kgl_sequence/kgl_sequence_codon.h"


namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coding Table Class
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Amino acid Translation tables
// Found at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
class AminoTranslationTable {

public:

  explicit AminoTranslationTable() : amino_table_rows_(*Tables::STANDARDTABLE) {}
  ~AminoTranslationTable() = default;

  [[nodiscard]] std::string TableName() const { return amino_table_rows_.table_name; }

  [[nodiscard]] std::string TableDescription() const { return amino_table_rows_.table_description; }

  // Get table by name.
  [[nodiscard]] bool setTranslationTable(const std::string& table_name);

  // Returns an amino acid for the codon. AminoAcid::Unknown if any bases are 'N'
  [[nodiscard]] AminoAcid::Alphabet getAmino(const Codon& codon);

  [[nodiscard]] bool isStopAmino(AminoAcid::Alphabet amino) const { return amino == AminoAcid::AMINO_STOP; }

  [[nodiscard]] bool isStartAmino(AminoAcid::Alphabet amino) const;

  [[nodiscard]] bool isStopCodon(const Codon& codon) const;

  [[nodiscard]] bool isStartCodon(const Codon& codon) const;

private:

  constexpr static size_t CONTAINS_BASE_N = 1000;

  TranslationTable amino_table_rows_;

  [[nodiscard]] size_t index(const Codon& Codon) const;

};



}   // end namespace

#endif //KGL_TABLE_H
