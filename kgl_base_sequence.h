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


class DNA5Sequence  {

public:

  using NucleotideType = typename NucleotideColumn_DNA5::NucleotideType ;
  using SequenceString = std::basic_string<NucleotideType>;

  explicit DNA5Sequence(SequenceString sequence) : base_sequence_(std::move(sequence)) {};
  DNA5Sequence() = delete;
  virtual ~DNA5Sequence() = default;

  NucleotideType operator[] (ContigOffset_t& offset) const { return base_sequence_[offset]; }
  ContigSize_t length() const { return base_sequence_.length(); }
  const NucleotideType* baseAddress(ContigOffset_t& offset) const { return &base_sequence_[offset]; }

  SequenceString getSequenceString() const { return base_sequence_; }


  static std::shared_ptr<DNA5Sequence> codingSequence(std::shared_ptr<const DNA5Sequence> base_sequence_ptr,
                                                      const SortedCDS& sorted_cds);

  std::shared_ptr<DNA5Sequence> codingSequence(const SortedCDS& sorted_cds) const {

    std::shared_ptr<const DNA5Sequence> seq_ptr(std::make_shared<const DNA5Sequence>(base_sequence_));
    return codingSequence(seq_ptr, sorted_cds);

  }

private:

  SequenceString base_sequence_;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_BASE_SEQUENCE_H
