//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SEQUENCE_BASE_H
#define KGL_SEQUENCE_BASE_H


#include <string>
#include "kgl_alphabet_base.h"
#include "kgl_genome_feature.h"
#include "kgl_sequence_virtual.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Base Sequence - A container for DNA/RNA sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using NucleotideType = Nucleotide_DNA5_t;
using SequenceString = std::basic_string<NucleotideType>;

class DNA5Sequence: public AlphabetSequence {

public:


  explicit DNA5Sequence(SequenceString sequence) : base_sequence_(std::move(sequence)) {};
  DNA5Sequence() = delete;
  ~DNA5Sequence() override = default;

  const std::string& getSequence() const override { return base_sequence_; }

  NucleotideType operator[] (ContigOffset_t& offset) const { return base_sequence_[offset]; }
  NucleotideType at(ContigOffset_t& offset) const { return base_sequence_[offset]; }

  ContigSize_t length() const { return base_sequence_.length(); }

  const NucleotideType* baseAddress(ContigOffset_t& offset) const { return &base_sequence_[offset]; }

  // Offset is the relative sequence offset.
  bool modifyBase(NucleotideType Nucleotide, ContigOffset_t sequence_offset);

  // Returns bool false if contig_offset is not within the coding sequences defined by sorted_cds.
  // If the contig_offset is in the coding sequence then a valid sequence_offset and the sequence length is returned.
  // The offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
  bool static offsetWithinSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                   ContigOffset_t contig_offset,
                                   ContigOffset_t& sequence_offset,
                                   ContigSize_t& sequence_length);

  // The contig_offset arg is used when then supplied sequence is a sub-sequence of the entire contig sequence
  // and is used to adjust calculated offsets when copying the CDS regions to the coding sequence.
  // For example, the supplied sub-sequence could be a Gene sequence, in this case contig_offset would be
  // the Gene offset within the contig. If the supplied sequence is the entire contig sequence then contig_offset
  // will be zero (the default argument value). The entire sequence defined by the sorted CDS is returned.
  static std::shared_ptr<DNA5Sequence> codingSequence(std::shared_ptr<const DNA5Sequence> contig_sequence_ptr,
                                                      std::shared_ptr<const CodingSequence> coding_seq_ptr) {

    return codingSubSequence(contig_sequence_ptr, coding_seq_ptr, 0, 0, 0);

  }

  std::shared_ptr<DNA5Sequence> codingSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr) const {

    return codingSubSequence(coding_seq_ptr, 0, 0);

  }

  // Same as the above, but returns a defined subsequence (generally a single/group of codons) of the coding sequence
  // Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the SortedCDS.
  static std::shared_ptr<DNA5Sequence> codingSubSequence(std::shared_ptr<const DNA5Sequence> base_sequence_ptr,
                                                         std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                         ContigOffset_t sub_sequence_offset, // base count; 0 == all.
                                                         ContigSize_t sub_sequence_length,  // number of bases; 0 == all
                                                         ContigOffset_t contig_offset = 0); // only if a subsequence.

  std::shared_ptr<DNA5Sequence> codingSubSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                  ContigOffset_t sub_sequence_offset,
                                                  ContigSize_t sub_sequence_length,
                                                  ContigOffset_t contig_offset = 0) const {

    std::shared_ptr<const DNA5Sequence> seq_ptr(std::make_shared<const DNA5Sequence>(base_sequence_));
    return codingSubSequence(seq_ptr, coding_seq_ptr, sub_sequence_offset, sub_sequence_length, contig_offset);

  }

private:

  SequenceString base_sequence_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_BASE_H
