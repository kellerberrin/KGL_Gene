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
// Base Sequence - Containers for DNA/RNA sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using NucleotideType = Nucleotide_DNA5_t;
using SequenceString = std::basic_string<NucleotideType>;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The actual sequence is contained in this base class. This also includes access routines.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DNA5Sequence: public AlphabetSequence {

public:


  explicit DNA5Sequence(SequenceString sequence) : base_sequence_(std::move(sequence)) {}
  DNA5Sequence() = delete;
  ~DNA5Sequence() override = default;


  NucleotideType operator[] (ContigOffset_t& offset) const { return base_sequence_[offset]; }
  NucleotideType at(ContigOffset_t& offset) const { return base_sequence_[offset]; }

  ContigSize_t length() const { return base_sequence_.length(); }

  // Offset is the relative sequence offset.
  bool modifyBase(NucleotideType Nucleotide, ContigOffset_t sequence_offset);


protected:

  SequenceString base_sequence_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a DNA5 coding sequence that has been generated using a CodingSequence (CDS) object.
// Only this object can be used to generate an amino acid sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DNA5SequenceCoding: public DNA5Sequence {

public:


  explicit DNA5SequenceCoding(SequenceString sequence) : DNA5Sequence(std::move(sequence)) {}
  DNA5SequenceCoding() = delete;
  ~DNA5SequenceCoding() override = default;

  std::string getSequenceAsString() const override { return static_cast<std::string>(base_sequence_); }

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that cannot be used to directly generate an Amino Acid sequence
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DNA5SequenceLinear: public DNA5Sequence {

public:


  explicit DNA5SequenceLinear(SequenceString sequence) : DNA5Sequence(std::move(sequence)) {}
  DNA5SequenceLinear() = delete;
  ~DNA5SequenceLinear() override = default;

  std::string getSequenceAsString() const override { return static_cast<std::string>(base_sequence_); }



  std::shared_ptr<DNA5SequenceCoding> codingSubSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                        ContigOffset_t sub_sequence_offset,
                                                        ContigSize_t sub_sequence_length,
                                                        ContigOffset_t contig_offset = 0) const {

    std::shared_ptr<const DNA5SequenceLinear> seq_ptr(std::make_shared<const DNA5SequenceLinear>(base_sequence_));
    return codingSubSequence(seq_ptr, coding_seq_ptr, sub_sequence_offset, sub_sequence_length, contig_offset);

  }

private:

  // Returns a defined subsequence (generally a single/group of codons) of the coding sequence
  // Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the SortedCDS.
  static std::shared_ptr<DNA5SequenceCoding> codingSubSequence(std::shared_ptr<const DNA5SequenceLinear> base_seq_ptr,
                                                               std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                               ContigOffset_t sub_sequence_offset, // base count
                                                               ContigSize_t sub_sequence_length,  // number of bases
                                                               ContigOffset_t contig_offset = 0); // if a subsequence.


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that cannot be used to directly generate an Amino Acid sequence
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DNA5SequenceContig: public DNA5SequenceLinear {

public:


  explicit DNA5SequenceContig(SequenceString sequence) : DNA5SequenceLinear(std::move(sequence)) {}
  DNA5SequenceContig() = delete;
  ~DNA5SequenceContig() override = default;


  // Returns bool false if contig_offset is not within the coding sequences defined by sorted_cds.
  // If the contig_offset is in the coding sequence then a valid sequence_offset and the sequence length is returned.
  // The offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
  bool offsetWithinSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                            ContigOffset_t contig_offset,
                            ContigOffset_t& sequence_offset,
                            ContigSize_t& sequence_length) const;

  // Returns the codon offset of offset within a coding, returns false if not within the coding sequence.
  bool codonOffset(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                   ContigOffset_t contig_offset,
                   ContigOffset_t& codon_offset,
                   ContigSize_t& base_in_codon) const;

  //The entire sequence defined by the sorted CDS is returned.
  std::shared_ptr<DNA5SequenceCoding> codingSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr) const {

    return codingSubSequence(coding_seq_ptr, 0, 0);

  }


private:


};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_BASE_H
