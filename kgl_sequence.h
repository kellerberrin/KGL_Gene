//
// Created by kellerberrin on 20/10/17.
//

#ifndef SAMFILE_KGL_SEQUENCE_H
#define SAMFILE_KGL_SEQUENCE_H

#include <string>
#include "kgl_nucleotide.h"
#include "kgl_genome_feature.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Convenience class to hold contig sequence strings.
template<typename T>
class BaseSequence {

public:

  explicit BaseSequence(std::basic_string<T> sequence) : base_sequence_(std::move(sequence)) {};
  BaseSequence() = delete;
  ~BaseSequence() = default;

  T operator[] (ContigOffset_t& offset) const { return base_sequence_[offset]; }
  ContigSize_t length() const { return base_sequence_.length(); }


  static BaseSequence codingSequence(std::shared_ptr<const BaseSequence> base_sequence_ptr,
                                     const SortedCDS& sorted_cds);

  BaseSequence codingSequence(const SortedCDS& sorted_cds) const {

    std::shared_ptr<const BaseSequence> seq_ptr(std::make_shared<const BaseSequence>(base_sequence_));
    return codingSequence(seq_ptr, sorted_cds);

  }


private:

  std::basic_string<T> base_sequence_;

};

using DNA5Sequence = BaseSequence<Nucleotide_DNA5_t>;

template<typename T>
BaseSequence<T> BaseSequence<T>::codingSequence(std::shared_ptr<const BaseSequence> base_sequence_ptr,
                                                const SortedCDS& sorted_cds) {

  // Check bounds.
  if (sorted_cds.rbegin()->second->sequence().end() >= base_sequence_ptr->length()) {

    ExecEnv::log().error("codingSequence(), CDS end offset: {} >= sequence size: {}",
                         sorted_cds.rbegin()->second->sequence().end(),
                         base_sequence_ptr->length());
    std::basic_string<T> null_seq;
    return BaseSequence<T>(null_seq);

  }

  // pre-allocate the string size.
  std::string::size_type seq_size = 0;
  for (auto cds : sorted_cds) {

    seq_size += cds.second->sequence().end() - cds.second->sequence().begin();

  }

  std::basic_string<T> coding_sequence;
  coding_sequence.reserve(seq_size);

  return BaseSequence<T>(coding_sequence);

}


// Given a contig sequence and an array of CDS, return the gene sequence string.
// If the strand is '-' then reverse and complement the sequence before returning.

/*
std::string sumof_filenames_1(std::vector<FileInfo> const& fv )
{
  std::string::size_type n = 0;
  for each ( FileInfo const& f in fv ) n +=  f.filename().size();

  std::string s;
  s.reserve( n );

  for( std::size_t ix = 0u; ix < fv.size(); ++ix) s += fv[ix].filename();

  return s;
}
*/


}   // namespace genome
}   // namespace kellerberrin

#endif //SAMFILE_KGL_SEQUENCE_H
