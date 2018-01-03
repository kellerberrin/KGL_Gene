//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SEQUENCE_BASE_H
#define KGL_SEQUENCE_BASE_H


#include <string>
#include <memory>
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

// A STRANDED DNA string that can be converted to an AMINO sequence.
class DNA5SequenceCoding: public AlphabetSequence<CodingDNA5> {

public:


  explicit DNA5SequenceCoding(StringCodingDNA5 sequence) : AlphabetSequence<CodingDNA5>(std::move(sequence)) {}
  explicit DNA5SequenceCoding(const DNA5SequenceCoding&) = default;
  DNA5SequenceCoding() = delete;
  ~DNA5SequenceCoding() override = default;

  // Offset is the relative sequence offset.
  bool modifyBase(ContigOffset_t sequence_offset, CodingDNA5::Alphabet Nucleotide);
  // Delete offset is relative to the begining of the sequence (0 is the first letter).
  bool deleteSubSequence(ContigOffset_t delete_offset, ContigSize_t delete_size);
  // Insert offset is relative to the begining of the sequence (0 is the first letter).
  bool insertSubSequence(ContigOffset_t insert_offset, const DNA5SequenceCoding& inserted_sequence);


private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that cannot be used to directly generate an Amino Acid sequence
// This sequence is NOT STRANDED.
// However, it can return the STRANDED sequence DNA5SequenceCoding using a CodingSequence (array of CDS).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A string of the standard 5 nucleotide DNA/RNA alphabet A, C, G, T/U, N
using StringDNA5 = AlphabetString<DNA5>;

// A sorted map of coding sequence offsets.
// The first (sort key) offset is the beginning of the EXON and the second (value) is the EXON end offset.
// The offsets use the half interval idiom [start, end).
// The offsets are initially contig offsets. However these can be modified to account for indel modifications.
// The map is generated by an adapter function exonOffsetAdapter() that uses a coding sequence as input.
using ExonOffsetMap = std::map<ContigOffset_t, ContigOffset_t>;

// An UNSTRANDED DNA string.
class DNA5SequenceLinear: public  AlphabetSequence<DNA5> {

public:


  explicit DNA5SequenceLinear(StringDNA5 sequence) :  AlphabetSequence<DNA5>(std::move(sequence)) {}
  explicit DNA5SequenceLinear() = default;
  explicit DNA5SequenceLinear(const DNA5SequenceLinear&) = default;
  ~DNA5SequenceLinear() override = default;

  // Returns a defined subsequence (generally a single/group of codons) of the coding sequence
  // Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the SortedCDS.
  std::shared_ptr<DNA5SequenceCoding> codingSubSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                        ContigOffset_t sub_sequence_offset,
                                                        ContigSize_t sub_sequence_length,
                                                        ContigOffset_t contig_offset = 0) const;
  // Offset is the relative sequence offset.
  bool modifyBase(ContigOffset_t base_offset, DNA5::Alphabet Nucleotide);
  // Delete offset is relative to the begining of the sequence (0 is the first letter).
  bool deleteSubSequence(ContigOffset_t delete_offset, ContigSize_t delete_size);
  // Insert offset is relative to the begining of the sequence (0 is the first letter).
  bool insertSubSequence(ContigOffset_t insert_offset, const DNA5SequenceLinear& inserted_sequence);

  std::string compareDNA5Sequences(const DNA5SequenceLinear& compare_seq) const { return compareSequences(compare_seq); }

private:

  // Returns a defined subsequence (generally a single/group of codons) of the coding sequence
  // Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the SortedCDS.
  static std::shared_ptr<DNA5SequenceCoding> codingSubSequence(std::shared_ptr<const DNA5SequenceLinear> base_seq_ptr,
                                                               const SortedCDS& sorted_cds,
                                                               StrandSense strand,
                                                               ContigOffset_t sub_sequence_offset, // base count
                                                               ContigSize_t sub_sequence_length,  // number of bases
                                                               ContigOffset_t contig_offset); // if a subsequence.

// A sorted map of coding sequence offsets.
// The first (sort key) offset is the beginning of the EXON and the second (value) is the EXON end offset.
// The offsets use the half interval idiom [start, end).
// The offsets are initially contig offsets. However these can be modified to account for indel offset modification.
// The map is generated by an adapter function exonOffsetAdapter() that uses a coding sequence as input.
  static bool exonOffsetAdapter(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                StrandSense& strand,
                                ExonOffsetMap& exon_offset_map);


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that implements a contig or chromosome.
// This sequence is NOT STRANDED.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UNSTRANDED
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

  // The entire sequence defined by the sorted CDS array is returned.
  std::shared_ptr<DNA5SequenceCoding> codingSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr) const {

    return codingSubSequence(coding_seq_ptr, 0, 0);

  }

  // Returns an UNSTRANDED region
  std::shared_ptr<DNA5SequenceLinear> unstrandedRegion(ContigOffset_t sub_sequence_offset, // offset
                                                       ContigSize_t sub_sequence_length) const { // if a subsequence.

    std::shared_ptr<DNA5SequenceLinear> sub_sequence(std::make_shared<DNA5SequenceLinear>());
    getSubsequence(sub_sequence_offset, sub_sequence_length, sub_sequence);
    return sub_sequence;

  }


private:


};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_BASE_H
