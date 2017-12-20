//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SEQUENCE_BASE_H
#define KGL_SEQUENCE_BASE_H


#include <string>
#include "kgl_alphabet_string.h"
#include "kgl_genome_feature.h"
#include "kgl_sequence_virtual.h"



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The DNA5 alphabet strings are defined here.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The actual sequence is contained in this base class. This also includes access routines.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a DNA5 coding sequence that has been generated using a CodingSequence (CDS) object.
// Only this object can be used to generate an amino acid sequence.
// This sequence is ALWAYS STRANDED.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A string of the standard 5 nucleotide DNA/RNA alphabet A, C, G, T/U, N
using StringCodingDNA5 = AlphabetString<CodingDNA5>;


class DNA5SequenceCoding: public AlphabetSequence<CodingDNA5> {

public:


  explicit DNA5SequenceCoding(StringCodingDNA5 sequence) : AlphabetSequence<CodingDNA5>(std::move(sequence)) {}
  DNA5SequenceCoding() = delete;
  ~DNA5SequenceCoding() override = default;

  // Offset is the relative sequence offset.
  bool modifyBase(CodingDNA5::Alphabet Nucleotide, ContigOffset_t sequence_offset);

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that cannot be used to directly generate an Amino Acid sequence
// This sequence is NOT STRANDED.
// However, it can return the STRANDED sequence DNA5SequenceCoding using a CodingSequence (array of CDS).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A string of the standard 5 nucleotide DNA/RNA alphabet A, C, G, T/U, N
using StringDNA5 = AlphabetString<DNA5>;


class DNA5SequenceLinear: public  AlphabetSequence<DNA5> {

public:


  explicit DNA5SequenceLinear(StringDNA5 sequence) :  AlphabetSequence<DNA5>(std::move(sequence)) {}
  DNA5SequenceLinear() = delete;
  ~DNA5SequenceLinear() override = default;

  std::shared_ptr<DNA5SequenceCoding> codingSubSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                        ContigOffset_t sub_sequence_offset,
                                                        ContigSize_t sub_sequence_length,
                                                        ContigOffset_t contig_offset = 0) const {

    std::shared_ptr<const DNA5SequenceLinear> seq_ptr(std::make_shared<const DNA5SequenceLinear>(alphabet_string_));
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
// A linear and contiguous DNA5 sequence that implements a contig or chromosome.
// This sequence is NOT STRANDED.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DNA5SequenceContig: public DNA5SequenceLinear {

public:


  explicit DNA5SequenceContig(StringDNA5 sequence) : DNA5SequenceLinear(std::move(sequence)) {}
  DNA5SequenceContig() = delete;
  ~DNA5SequenceContig() override = default;


  // Returns bool false if contig_offset is not within the coding sequences defined by sorted_cds.
  // If the contig_offset is in the coding sequence then a valid sequence_offset and the sequence length is returned.
  // The offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
  bool offsetWithinSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                            ContigOffset_t contig_offset,
                            ContigOffset_t& sequence_offset,
                            ContigSize_t& sequence_length) const;

  // Returns the codon offset of the contig offset within a coding sequence, returns false if not within the coding sequence.
  bool codonOffset(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                   ContigOffset_t contig_offset,
                   ContigOffset_t& codon_offset,
                   ContigSize_t& base_in_codon) const;

  //The entire sequence defined by the sorted CDS array is returned.
  std::shared_ptr<DNA5SequenceCoding> codingSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr) const {

    return codingSubSequence(coding_seq_ptr, 0, 0);

  }


private:


};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_BASE_H
