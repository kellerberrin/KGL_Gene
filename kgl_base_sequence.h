//
// Created by kellerberrin on 22/10/17.
//

#ifndef KGL_BASE_SEQUENCE_H
#define KGL_BASE_SEQUENCE_H


#include <string>
#include "kgl_nucleotide.h"
#include "kgl_genome_feature.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Base Sequence - A container for DNA/RNA sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T>
class BaseSequence  {

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
  std::string::size_type calculated_seq_size = 0;
  for (auto cds : sorted_cds) {

    calculated_seq_size += cds.second->sequence().end() - cds.second->sequence().begin();

  }

  SequenceString coding_sequence;
  coding_sequence.reserve(calculated_seq_size + 1); // Just to make sure.

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
      break;

    case StrandSense::REVERSE: {

      // Insert in reverse complement order.
      typename SequenceString::const_reverse_iterator rbegin;
      typename SequenceString::const_reverse_iterator rend;
      auto complement_base = [](typename BaseSequence::NucleotideType base) { return T::complementNucleotide(base); };
      for (auto rit = sorted_cds.rbegin(); rit != sorted_cds.rend(); ++rit) {

        rbegin = base_sequence_ptr->base_sequence_.rbegin()
                 + (base_sequence_ptr->length() - rit->second->sequence().end());
        rend = base_sequence_ptr->base_sequence_.rbegin()
               + (base_sequence_ptr->length() - rit->second->sequence().begin());
        std::transform( rbegin, rend, std::back_inserter(coding_sequence), complement_base);

      }


    }
      break;

      // Check the sequence size.

      if (coding_sequence.length() != calculated_seq_size) {

        ExecEnv::log().error("Coding sequence length: {} NOT EQUAL to calculated sequence: {}",
                             coding_sequence.length(),
                             calculated_seq_size);

      }

  }

  return std::shared_ptr<BaseSequence<T>>(std::make_shared<BaseSequence<T>>(BaseSequence(coding_sequence)));

}


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_BASE_SEQUENCE_H
