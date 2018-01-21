//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SEQUENCE_AMINO_H
#define KGL_SEQUENCE_AMINO_H


#include <string>
#include "kgl_sequence_base.h"
#include "kgl_table.h"
#include "kgl_sequence_virtual.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Amino alphabet strings are defined here.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using StringAminoAcid = AlphabetString<AminoAcid>;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino Sequence - A container for Amino Acid (protein) sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class AminoSequence: public AlphabetSequence<AminoAcid> {

public:

  explicit AminoSequence(StringAminoAcid sequence) : AlphabetSequence<AminoAcid>(std::move(sequence)) {};
  AminoSequence() = delete;
  ~AminoSequence() override = default;

  std::string compareAminoSequences(std::shared_ptr<const AminoSequence> compare_seq, CompareScore_t& score) const {

    return compareSequencesAmino(*compare_seq, score);

  }

  static std::string multipleCompare(const std::vector<std::shared_ptr<const AminoSequence>>& compare_seq_vec);


  bool removeTrailingStop();  // Remove the stop codon (if present).

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TranslateToAmino - Convert DNA/RNA base sequences to Amino acid sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TranslateToAmino {

public:


  TranslateToAmino() : table_ptr_(std::make_shared<AminoTranslationTable>())  {}
  ~TranslateToAmino() = default;

  std::string translationTableName() const { return table_ptr_->TableName(); }

  std::string translationTableDescription() const { return table_ptr_->TableDescription(); }

  bool settranslationTable(const std::string& table_name) { return table_ptr_->setTranslationTable(table_name); }

  Codon firstCodon(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const {

    return Codon(sequence_ptr, 0);

  }

  bool checkStartCodon(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const {

    return table_ptr_->isStartCodon(firstCodon(sequence_ptr));

  }

  bool checkStartCodon(std::shared_ptr<const AminoSequence> sequence_ptr) const;
  bool checkStopCodon(std::shared_ptr<const AminoSequence> sequence_ptr) const;
  size_t checkNonsenseMutation(std::shared_ptr<const AminoSequence> sequence_ptr) const;

  std::shared_ptr<AminoSequence> getAminoSequence(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const;

  std::shared_ptr<AminoSequence> getAminoSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                  std::shared_ptr<const DNA5SequenceContig> contig_sequence_ptr) const;

  Codon lastCodon(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const {

    return Codon(sequence_ptr, Codon::codonLength(sequence_ptr) - 1);

  }

  bool checkStopCodon(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const {

    return table_ptr_->isStopCodon(lastCodon(sequence_ptr));

  }

  size_t checkNonsenseMutation(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const;

  AminoAcid::Alphabet getAmino(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr, ContigSize_t codon_index) const;
  AminoAcid::Alphabet getAmino(const Codon& codon) const;


private:

  std::shared_ptr<AminoTranslationTable> table_ptr_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_AMINO_H
