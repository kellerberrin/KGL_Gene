//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SEQUENCE_AMINO_H
#define KGL_SEQUENCE_AMINO_H

#include "kgl_sequence_base_view.h"
#include "kgl_sequence_base.h"
#include "kgl_table.h"

#include <memory>

namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Amino alphabet strings are defined here.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


using StringAminoAcid = AlphabetString<AminoAcid>;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino Sequence - A container for Amino Acid (protein) sequences.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AminoSequence: public AlphabetSequence<AminoAcid> {

public:

  AminoSequence(AminoSequence&& sequence) noexcept : AlphabetSequence<AminoAcid>(std::move(sequence)) {}
  explicit AminoSequence(const AminoSequenceView& sequence_view) : AlphabetSequence<AminoAcid>(sequence_view.getSequence()) {}
  explicit AminoSequence(StringAminoAcid&& sequence_string) noexcept : AlphabetSequence<AminoAcid>(std::move(sequence_string)) {}
  AminoSequence(const AminoSequence& sequence) = delete; // For Performance reasons, don't allow copy constructors.
  ~AminoSequence() override = default;

  // For Performance reasons, don't allow naive assignments.
  AminoSequence& operator=(const AminoSequence&) = delete;
  // Only allow move assignments
  AminoSequence& operator=(AminoSequence&& moved) noexcept {

    AlphabetSequence<AminoAcid>::operator=(std::move(moved));
    return *this;

  }

  ///Return a view.
  [[nodiscard]] AminoSequenceView getView() const { return AminoSequenceView(*this); }

private:


};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TranslateToAmino - Convert DNA/RNA base sequences to Amino acid sequences.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TranslateToAmino {

public:


  TranslateToAmino() = default;
  ~TranslateToAmino() = default;

  /// ReturnType the name of the translation table.
  [[nodiscard]] std::string translationTableName() const { return table_.TableName(); }

  /// Set the translation table by name.
  [[nodiscard]] bool settranslationTable(const std::string& table_name) { return table_.setTranslationTable(table_name); }


  [[nodiscard]] bool checkStartCodon(const AminoSequence& amino_sequence) const;
  [[nodiscard]] bool checkStopCodon(const AminoSequence& amino_sequence) const;
  /// Find the size of the sequence including the first stop codon.
  /// If not found, bool is false and the sequence length is returned.
  [[nodiscard]] std::pair<size_t, bool> firstStopSequenceSize(const AminoSequence& amino_sequence) const;

  /// ReturnType the amino acid sequence for the coding sequence.
  [[nodiscard]] AminoSequence getAminoSequence(const DNA5SequenceCoding& coding_sequence) const;
  /// ReturnType the amino acid for the codon. AminoAcid::Unknown if any bases are 'N'
  [[nodiscard]] AminoAcid::Alphabet getAmino(const Codon& codon) const;


private:

  // The translation table is never null and never shared: a value member removes
  // the shared_ptr indirection on the per-codon translation hot path.
  AminoTranslationTable table_;

};


}   // end namespace


#endif //KGL_SEQUENCE_AMINO_H