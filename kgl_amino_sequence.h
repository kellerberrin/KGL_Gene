//
// Created by kellerberrin on 29/10/17.
//

#ifndef KGL_AMINO_SEQUENCE_H
#define KGL_AMINO_SEQUENCE_H


#include <string>
#include "kgl_base_sequence.h"
#include "kgl_table.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino Sequence - A container for Amino Acid (protein) sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AminoSequence {

public:

  using AminoType = typename AminoAcidTypes::AminoType;
  using ProteinString = std::basic_string<AminoType>;

  explicit AminoSequence(ProteinString sequence) : amino_sequence_(std::move(sequence)) {};
  AminoSequence() = delete;
  virtual ~AminoSequence() = default;

  AminoType operator[] (ContigOffset_t& offset) const { return amino_sequence_[offset]; }
  ContigSize_t length() const { return amino_sequence_.length(); }


private:

  ProteinString amino_sequence_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coding Sequence - Convert DNA/RNA base sequences to Amino acid sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CodingSequenceDNA5 {

public:


  CodingSequenceDNA5() : table_ptr_(std::make_shared<AminoTranslationTable>())  {}
  ~CodingSequenceDNA5() = default;

  std::string translationTableName() const { return table_ptr_->TableName(); }

  bool settranslationTable(size_t table) { return table_ptr_->setTranslationTable(table); }

  AminoAcidTypes::Codon firstCodon(std::shared_ptr<DNA5Sequence> sequence_ptr) const {

    return getCodon(sequence_ptr, 0);

  }

  bool checkStartCodon(std::shared_ptr<DNA5Sequence> sequence_ptr) const {

    return table_ptr_->isStartCodon(firstCodon(sequence_ptr));

  }

  std::shared_ptr<AminoSequence> getAminoSequence(std::shared_ptr<DNA5Sequence> sequence_ptr) const;

  AminoAcidTypes::Codon lastCodon(std::shared_ptr<DNA5Sequence> sequence_ptr) const {

    return getCodon(sequence_ptr, codonLength(sequence_ptr) - 1);

  }

  bool checkStopCodon(std::shared_ptr<DNA5Sequence> sequence_ptr) const {

    return table_ptr_->isStopCodon(lastCodon(sequence_ptr));

  }

  size_t checkNonsenseMutation(std::shared_ptr<DNA5Sequence> sequence_ptr) const;

  size_t codonLength(std::shared_ptr<DNA5Sequence> sequence_ptr) const { return sequence_ptr->length() / 3; }

  AminoAcidTypes::Codon getCodon(std::shared_ptr<DNA5Sequence> sequence_ptr, size_t index) const;

private:

  std::shared_ptr<AminoTranslationTable> table_ptr_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //READSAMFILE_KGL_AMINO_SEQUENCE_H
