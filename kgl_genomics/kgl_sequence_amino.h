//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SEQUENCE_AMINO_H
#define KGL_SEQUENCE_AMINO_H


#include <string>
#include "kgl_sequence_base.h"
#include "kgl_table.h"
#include "kgl_sequence_virtual.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Amino alphabet strings are defined here.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using StringAminoAcid = AlphabetString<AminoAcid>;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino Sequence - A container for Amino Acid (protein) sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class AminoSequence: public AlphabetSequence<AminoAcid> {

public:

  AminoSequence(AminoSequence&& sequence) noexcept : AlphabetSequence<AminoAcid>(std::move(sequence)) {};
  explicit AminoSequence(StringAminoAcid&& sequence_string) noexcept : AlphabetSequence<AminoAcid>(std::move(sequence_string)) {};
  AminoSequence(AminoSequence& sequence) = delete; // For Performance reasons, don't allow copy constructors.
  AminoSequence() = default;  // Allow the creation of empty sequences
  ~AminoSequence() override = default;

  // For Performance reasons, don't allow naive assignments.
  AminoSequence& operator=(const AminoSequence&) = delete;
  // Only allow move assignments
  AminoSequence& operator=(AminoSequence&& moved) noexcept {

    AlphabetSequence<AminoAcid>::operator=(std::move(moved));
    return *this;

  }

  [[nodiscard]] bool removeTrailingStop();  // Remove the stop codon (if present).

  // Returns a sub-sequence
  [[nodiscard]] AminoSequence subSequence( ContigOffset_t sub_sequence_offset, ContigSize_t sub_sequence_length) const;

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TranslateToAmino - Convert DNA/RNA base sequences to Amino acid sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TranslateToAmino {

public:


  TranslateToAmino() : table_ptr_(std::make_shared<AminoTranslationTable>())  {}
  ~TranslateToAmino() = default;

  [[nodiscard]] std::string translationTableName() const { return table_ptr_->TableName(); }

  [[nodiscard]] std::string translationTableDescription() const { return table_ptr_->TableDescription(); }

  [[nodiscard]] bool settranslationTable(const std::string& table_name) { return table_ptr_->setTranslationTable(table_name); }

  [[nodiscard]] Codon firstCodon(const DNA5SequenceCoding& coding_sequence) const {

    return Codon(coding_sequence, 0);

  }

  [[nodiscard]] bool checkStartCodon(const DNA5SequenceCoding& coding_sequence) const {

    return table_ptr_->isStartCodon(firstCodon(coding_sequence));

  }

  [[nodiscard]] bool checkStartCodon(const AminoSequence& amino_sequence) const;
  [[nodiscard]] bool checkStopCodon(const AminoSequence& amino_sequence) const;
  [[nodiscard]] size_t checkNonsenseMutation(const AminoSequence& amino_sequence) const;

  [[nodiscard]] AminoSequence getAminoSequence(const DNA5SequenceCoding& coding_sequence) const;

  [[nodiscard]] AminoSequence getAminoSequence( const std::shared_ptr<const CodingSequence>& coding_seq_ptr,
                                                const DNA5SequenceContig& contig_sequence) const;

  [[nodiscard]] Codon lastCodon(const DNA5SequenceCoding& coding_sequence) const {

    return Codon(coding_sequence, Codon::codonLength(coding_sequence) - 1);

  }

  [[nodiscard]] bool checkStopCodon(const DNA5SequenceCoding& coding_sequence) const {

    return table_ptr_->isStopCodon(lastCodon(coding_sequence));

  }

  [[nodiscard]] size_t checkNonsenseMutation(const DNA5SequenceCoding& coding_sequence) const;

  [[nodiscard]] AminoAcid::Alphabet getAmino(const DNA5SequenceCoding& coding_sequence, ContigSize_t codon_index) const;
  [[nodiscard]] AminoAcid::Alphabet getAmino(const Codon& codon) const;


private:

  std::shared_ptr<AminoTranslationTable> table_ptr_;

};


}   // end namespace


#endif //KGL_SEQUENCE_AMINO_H
