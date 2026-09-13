//
// Created by kellerberrin on 22/10/17.
//

#ifndef KGL_TABLE_H
#define KGL_TABLE_H


#include "kgl_table_registry.h"
#include "kgl_sequence_codon.h"


namespace kellerberrin::genome {   //  organization level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coding Table Class
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Amino acid Translation tables
// Found at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
class AminoTranslationTable {

public:

  AminoTranslationTable() : amino_table_rows_(*Tables::STANDARDTABLE) {}
  ~AminoTranslationTable() = default;

  /// ReturnType the translation table name.
  [[nodiscard]] std::string TableName() const { return amino_table_rows_.table_name; }

  /// Get table by name.
  [[nodiscard]] bool setTranslationTable(const std::string& table_name);

  /// Returns an amino acid for the codon. AminoAcid::Unknown if any bases are 'N'
  [[nodiscard]] AminoAcid::Alphabet getAmino(const Codon& codon) const;

  /// Returns true if the amino acid is a stop amino acid.
  [[nodiscard]] bool isStopAmino(AminoAcid::Alphabet amino) const { return amino == AminoAcid::AMINO_STOP; }

  /// Returns true if any codon in the table that codes the amino acid is a start codon.
  [[nodiscard]] bool isStartAmino(AminoAcid::Alphabet amino) const;

  /// Returns true if the codon is a stop codon.
  [[nodiscard]] bool isStopCodon(const Codon& codon) const;

  /// Returns true if the codon is a start codon.
  [[nodiscard]] bool isStartCodon(const Codon& codon) const;

private:

  constexpr static size_t CONTAINS_BASE_N = 1000;

  TranslationTable amino_table_rows_;

  /// ReturnType the table index for the codon or CONTAINS_BASE_N if the codon contains an 'N' base.
  [[nodiscard]] size_t index(const Codon& codon) const;

};


}   // end namespace

#endif //KGL_TABLE_H