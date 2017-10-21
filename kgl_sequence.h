//
// Created by kellerberrin on 20/10/17.
//

#ifndef SAMFILE_KGL_SEQUENCE_H
#define SAMFILE_KGL_SEQUENCE_H

#include <string>
#include "kgl_nucleotide.h"
#include "kgl_amino.h"
#include "kgl_genome_feature.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// Convenience class to hold contig sequence strings.
template<typename T>
class BaseSequence {

public:

  using NucleotideType = typename T::NucleotideType;
  using SequenceString = std::basic_string<NucleotideType>;

  explicit BaseSequence(SequenceString sequence) : base_sequence_(std::move(sequence)) {};
  BaseSequence() = delete;
  virtual ~BaseSequence() = default;

  NucleotideType operator[] (ContigOffset_t& offset) const { return base_sequence_[offset]; }
  ContigSize_t length() const { return base_sequence_.length(); }
  const NucleotideType* baseAddress(ContigOffset_t& offset) const { return &base_sequence_[offset]; }


  static std::shared_ptr<BaseSequence> codingSequence(std::shared_ptr<const BaseSequence> base_sequence_ptr,
                                                      const SortedCDS& sorted_cds);

  std::shared_ptr<BaseSequence> codingSequence(const SortedCDS& sorted_cds) const {

    std::shared_ptr<const BaseSequence> seq_ptr(std::make_shared<const BaseSequence>(base_sequence_));
    return codingSequence(seq_ptr, sorted_cds);

  }


private:

  SequenceString base_sequence_;

};


template<class T, class AminoTable>
class CodingSequence {

public:


  explicit CodingSequence(std::shared_ptr<BaseSequence<T>> sequence_ptr) : sequence_ptr_(sequence_ptr) {};
  CodingSequence() = delete;
  ~CodingSequence() = default;

  inline AminoAcidTypes::Codon firstCodon() const { return getCodon(0); }
  bool checkStartCodon() const { return AminoTable::isStartCodon(firstCodon()); }
  inline AminoAcidTypes::Codon lastCodon() const { return getCodon(codonLength() - 1); }
  bool checkStopCodon() const { return AminoTable::isStopCodon(lastCodon()); }
  size_t checkNonsenseMutation() const {

    for (size_t index = 0; index < codonLength() - 1; ++index) {

      if (AminoTable::isStopCodon(getCodon(index))) return index;

    }

    return 0;

  }
  inline size_t codonLength() const { return sequence_ptr_->length() / 3; }
  inline AminoAcidTypes::Codon getCodon(size_t index) const {

    if (index >= codonLength()) {

      ExecEnv::log().error("Invalid codon specified index:{}, for coding sequence length:{} (first codon returned)",
                           index, sequence_ptr_->length());
      index = 0;
    }

    AminoAcidTypes::Codon codon;
    index = index * 3;
    codon.bases = sequence_ptr_->baseAddress(index);
    return codon;

  }

private:

  std::shared_ptr<BaseSequence<T>> sequence_ptr_;

};


using DNA5Sequence = BaseSequence<NucleotideColumn_DNA5>;
using StandardCodingSequence = CodingSequence<NucleotideColumn_DNA5, StandardAminoTranslationTable>;


template<typename T>
std::shared_ptr<BaseSequence<T>> BaseSequence<T>::codingSequence(std::shared_ptr<const BaseSequence> base_sequence_ptr,
                                                const SortedCDS& sorted_cds) {

  // If no cds then return null string.
  if (sorted_cds.empty()) {

    SequenceString null_seq;
    return std::shared_ptr<BaseSequence<T>>(std::make_shared<BaseSequence<T>>(BaseSequence(null_seq)));

  }

  // Check bounds.
  if (sorted_cds.rbegin()->second->sequence().end() >= base_sequence_ptr->length()) {

    ExecEnv::log().error("codingSequence(), CDS end offset: {} >= target sequence size: {}",
                         sorted_cds.rbegin()->second->sequence().end(),
                         base_sequence_ptr->length());
    SequenceString null_seq;
    return std::shared_ptr<BaseSequence<T>>(std::make_shared<BaseSequence<T>>(BaseSequence(null_seq)));

  }

  // For efficiency, pre-allocate storage for the sequence string.
  std::string::size_type seq_size = 0;
  for (auto cds : sorted_cds) {

    seq_size += cds.second->sequence().end() - cds.second->sequence().begin();

  }

  SequenceString coding_sequence;
  coding_sequence.reserve(seq_size + 1); // Just to make sure.

  // Get the strand and copy or reverse copy the base complement.
  switch(sorted_cds.begin()->second->sequence().strand()) {

    case StrandSense::UNKNOWN:
      ExecEnv::log().warn("codingSequence(); CDS: {} offset: {} has 'UNKNOWN' ('.') strand assuming 'FORWARD' ('+')",
                          sorted_cds.begin()->second->id(),
                          sorted_cds.begin()->second->sequence().begin());
    case StrandSense::FORWARD: {

      typename SequenceString::const_iterator begin;
      typename SequenceString::const_iterator end;
      for (auto cds : sorted_cds) {

        begin = base_sequence_ptr->base_sequence_.begin() + cds.second->sequence().begin();
        end = base_sequence_ptr->base_sequence_.begin() + cds.second->sequence().end();
        std::copy( begin, end, std::back_inserter(coding_sequence));

      }

    }

    case StrandSense::REVERSE: {

      // Insert in reverse complement order.
      typename SequenceString::const_reverse_iterator rbegin;
      typename SequenceString::const_reverse_iterator rend;
      auto complement_base = [](typename BaseSequence::NucleotideType base) { return T::complementNucleotide(base); };
      for (auto cds : sorted_cds) {

        rbegin = base_sequence_ptr->base_sequence_.rbegin()
                 + (base_sequence_ptr->length() - cds.second->sequence().end());
        rend = base_sequence_ptr->base_sequence_.rbegin()
               + (base_sequence_ptr->length() - cds.second->sequence().begin());
        std::transform( rbegin, rend, std::back_inserter(coding_sequence), complement_base);

      }

    }

  }

  return std::shared_ptr<BaseSequence<T>>(std::make_shared<BaseSequence<T>>(BaseSequence(coding_sequence)));

}


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_SEQUENCE_H
