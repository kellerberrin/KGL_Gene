//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SEQUENCE_AMINO_H
#define KGL_SEQUENCE_AMINO_H


#include <string>
#include "kgl_sequence_base.h"
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

  const ProteinString& getProteinString() const { return amino_sequence_; }

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

  std::shared_ptr<AminoSequence> getAminoSequence(const SortedCDS& sorted_cds,
                                                  std::shared_ptr<const DNA5Sequence> contig_sequence_ptr) const;

  AminoAcidTypes::Codon lastCodon(std::shared_ptr<DNA5Sequence> sequence_ptr) const {

    return getCodon(sequence_ptr, codonLength(sequence_ptr) - 1);

  }

  bool checkStopCodon(std::shared_ptr<DNA5Sequence> sequence_ptr) const {

    return table_ptr_->isStopCodon(lastCodon(sequence_ptr));

  }

  size_t checkNonsenseMutation(std::shared_ptr<DNA5Sequence> sequence_ptr) const;

  static ContigOffset_t codonLength(std::shared_ptr<DNA5Sequence> sequence_ptr) {

    return static_cast<ContigOffset_t>(sequence_ptr->length() / 3);

  }

  static AminoAcidTypes::Codon getCodon(std::shared_ptr<DNA5Sequence> sequence_ptr, ContigOffset_t codon_index);

  AminoAcidTypes::AminoType getAmino(std::shared_ptr<DNA5Sequence> sequence_ptr, ContigSize_t codon_index) const;

  static bool codonOffset(const SortedCDS& sorted_cds,
                          ContigOffset_t contig_offset,
                          ContigOffset_t& codon_offset,
                          ContigSize_t& base_in_codon);

  // Returns the amino mutation of an SNP in a coding sequence.
  bool SNPMutation(const SortedCDS& sorted_cds,
                   const std::shared_ptr<const DNA5Sequence>& contig_sequence_ptr,
                   ContigOffset_t contig_offset,
                   typename NucleotideColumn_DNA5::NucleotideType reference_base,
                   typename NucleotideColumn_DNA5::NucleotideType mutant_base,
                   ContigOffset_t& codon_offset,
                   typename AminoAcidTypes::AminoType& reference_amino,
                   typename AminoAcidTypes::AminoType& mutant_amino) const;

private:

  std::shared_ptr<AminoTranslationTable> table_ptr_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_AMINO_H
