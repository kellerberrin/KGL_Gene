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
// The Amino alphabet strings are defined here.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using StringAminoAcid = AlphabetString<AminoAcid>;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino Sequence - A container for Amino Acid (protein) sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class AminoSequence: public AlphabetSequence {

public:

  explicit AminoSequence(StringAminoAcid sequence) : amino_sequence_(std::move(sequence)) {};
  AminoSequence() = delete;
  ~AminoSequence() override = default;

  std::string getSequenceAsString() const override { return amino_sequence_.str(); }

  AminoAcid::Alphabet operator[] (ContigOffset_t& offset) const { return amino_sequence_[offset]; }
  ContigSize_t length() const { return amino_sequence_.length(); }

  const std::string compareAminoSequences(const std::shared_ptr<const AminoSequence>& amino_seq_ptr) const {

    return compareAminoSequences(amino_seq_ptr->amino_sequence_, amino_sequence_);

  }

  const std::string emphasizeAminoSequence(const std::vector<ContigOffset_t>& emphasize_offsets) const {

    return emphasizeProteinString(amino_sequence_, emphasize_offsets);

  }

  bool removeTrailingStop();  // Remove the stop codon (if present).

private:

  StringAminoAcid amino_sequence_;

  static std::string emphasizeProteinString(const StringAminoAcid& amino_string,
                                            const std::vector<ContigOffset_t>& emphasize_offsets);

  static std::string compareAminoSequences(const StringAminoAcid& compare_amino,
                                           const StringAminoAcid& reference_amino);

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TranslateToAmino - Convert DNA/RNA base sequences to Amino acid sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TranslateToAmino {

public:


  TranslateToAmino() : table_ptr_(std::make_shared<AminoTranslationTable>())  {}
  ~TranslateToAmino() = default;

  std::string translationTableName() const { return table_ptr_->TableName(); }

  bool settranslationTable(size_t table) { return table_ptr_->setTranslationTable(table); }

  Codon firstCodon(std::shared_ptr<DNA5SequenceCoding> sequence_ptr) const {

    return Codon(sequence_ptr, 0);

  }

  bool checkStartCodon(std::shared_ptr<DNA5SequenceCoding> sequence_ptr) const {

    return table_ptr_->isStartCodon(firstCodon(sequence_ptr));

  }

  std::shared_ptr<AminoSequence> getAminoSequence(std::shared_ptr<DNA5SequenceCoding> sequence_ptr) const;

  std::shared_ptr<AminoSequence> getAminoSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                  std::shared_ptr<const DNA5SequenceContig> contig_sequence_ptr) const;

  Codon lastCodon(std::shared_ptr<DNA5SequenceCoding> sequence_ptr) const {

    return Codon(sequence_ptr, Codon::codonLength(sequence_ptr) - 1);

  }

  bool checkStopCodon(std::shared_ptr<DNA5SequenceCoding> sequence_ptr) const {

    return table_ptr_->isStopCodon(lastCodon(sequence_ptr));

  }

  size_t checkNonsenseMutation(std::shared_ptr<DNA5SequenceCoding> sequence_ptr) const;

  AminoAcid::Alphabet getAmino(std::shared_ptr<DNA5SequenceCoding> sequence_ptr, ContigSize_t codon_index) const;

  // Returns the amino mutation of an SNP in a coding sequence.
  bool SNPMutation(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                   const std::shared_ptr<const DNA5SequenceContig>& contig_sequence_ptr,
                   ContigOffset_t contig_offset,
                   DNA5::Alphabet reference_base,
                   DNA5::Alphabet mutant_base,
                   ContigOffset_t& codon_offset,
                   AminoAcid::Alphabet& reference_amino,
                   AminoAcid::Alphabet& mutant_amino) const;

private:

  std::shared_ptr<AminoTranslationTable> table_ptr_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_AMINO_H
