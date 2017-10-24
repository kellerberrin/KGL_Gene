//
// Created by kellerberrin on 20/10/17.
//

#ifndef SAMFILE_KGL_SEQUENCE_H
#define SAMFILE_KGL_SEQUENCE_H

#include <string>
#include "kgl_base_sequence.h"
#include "kgl_table.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino Sequence - A container for Amino Acid (protein) sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
class AminoSequence {

public:

  using AminoType = typename T::AminoType;
  using ProteinString = std::basic_string<AminoType>;

  explicit AminoSequence(ProteinString sequence) : amino_sequence_(std::move(sequence)) {};
  AminoSequence() = delete;
  virtual ~AminoSequence() = default;

  AminoType operator[] (ContigOffset_t& offset) const { return amino_sequence_[offset]; }
  ContigSize_t length() const { return amino_sequence_.length(); }
  const AminoType* baseAddress(ContigOffset_t& offset) const { return &amino_sequence_[offset]; }


private:

  ProteinString amino_sequence_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino Sequence Specialization for the Standard Amino Acid types.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using StandardAminoSequence = AminoSequence<AminoAcidTypes>;

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

  // The NCBI Amino Acid translation tables (1-31).
  bool setTranslationTable(size_t table) {

    bool result = false;

    if(table > Tables::TABLE_ARRAY_SIZE) {

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

  AminoAcidTypes::AminoType getAmino(const AminoAcidTypes::Codon& Codon) {

    return amino_table_rows_.amino_table[index(Codon)].amino_acid;

  }

  bool isStopCodon(const AminoAcidTypes::Codon& Codon) {

    return (amino_table_rows_.amino_table[index(Codon)].start == AminoAcidTypes::STOP_CODON);

  }

  bool isStartCodon(const AminoAcidTypes::Codon& Codon) {

    return amino_table_rows_.amino_table[index(Codon)].start == AminoAcidTypes::START_CODON;

  }

private:

  TranslationTable amino_table_rows_;

  inline size_t index(const AminoAcidTypes::Codon& Codon) {

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

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coding Sequence - An assembly of DNA/RNA sequences from CDS features.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
class CodingSequence {

public:


  CodingSequence() : table_ptr_(std::make_shared<AminoTranslationTable>())  {}
  ~CodingSequence() = default;

  std::string translationTableName() const { return table_ptr_->TableName(); }

  bool settranslationTable(size_t table) { return table_ptr_->setTranslationTable(table); }

  AminoAcidTypes::Codon firstCodon(std::shared_ptr<BaseSequence<T>> sequence_ptr) const {

    return getCodon(sequence_ptr, 0);

  }

  bool checkStartCodon(std::shared_ptr<BaseSequence<T>> sequence_ptr) const {

    return table_ptr_->isStartCodon(firstCodon(sequence_ptr));

  }

  AminoAcidTypes::Codon lastCodon(std::shared_ptr<BaseSequence<T>> sequence_ptr) const
  {

    return getCodon(sequence_ptr, codonLength(sequence_ptr) - 1);

  }

  bool checkStopCodon(std::shared_ptr<BaseSequence<T>> sequence_ptr) const {

    return table_ptr_->isStopCodon(lastCodon(sequence_ptr));

  }

  size_t checkNonsenseMutation(std::shared_ptr<BaseSequence<T>> sequence_ptr) const {

    for (size_t index = 0; index < codonLength(sequence_ptr) - 1; ++index) {

      if (table_ptr_->isStopCodon(getCodon(sequence_ptr,index))) return index;

    }

    return 0;

  }

  inline size_t codonLength(std::shared_ptr<BaseSequence<T>> sequence_ptr) const { return sequence_ptr->length() / 3; }

  inline AminoAcidTypes::Codon getCodon(std::shared_ptr<BaseSequence<T>> sequence_ptr, size_t index) const {

    if (index >= codonLength(sequence_ptr)) {

      ExecEnv::log().error("Invalid codon specified index:{}, for coding sequence length:{} (first codon returned)",
                           index, sequence_ptr->length());
      index = 0;
    }

    AminoAcidTypes::Codon codon;
    index = index * 3;
    codon.bases = sequence_ptr->baseAddress(index);
    return codon;

  }

private:

  std::shared_ptr<AminoTranslationTable> table_ptr_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coding Sequence Specialization for DNA
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using CodingSequenceDNA5 = CodingSequence<NucleotideColumn_DNA5>;


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_SEQUENCE_H
