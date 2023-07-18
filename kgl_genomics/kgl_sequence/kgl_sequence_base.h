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



namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The DNA5 alphabet strings are defined here.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The actual sequence is contained in the base class. This also includes access routines.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a DNA5 coding sequence that has been generated using a TranscriptionSequence (CDS) object.
// Or upconverted from a DNA5SequenceLinear sequence (see function codingSequence(StrandSense strand) below).
// Only this object can be used to generate an amino acid sequence.
// This sequence is ALWAYS STRANDED.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A string of the standard 5 nucleotide DNA/RNA alphabet A, C, G, T/U, N
using StringCodingDNA5 = AlphabetString<CodingDNA5>;

// A STRANDED DNA string that can be converted to an AMINO sequence.
class DNA5SequenceCoding: public AlphabetSequence<CodingDNA5> {

public:

  DNA5SequenceCoding(DNA5SequenceCoding&& sequence) noexcept : AlphabetSequence<CodingDNA5>(std::move(sequence)), strand_(sequence.strand_) {}
  explicit DNA5SequenceCoding(StringCodingDNA5&& sequence_string, StrandSense strand) noexcept : AlphabetSequence<CodingDNA5>(std::move(sequence_string)), strand_(strand) {}
  DNA5SequenceCoding() : strand_(StrandSense::FORWARD) {} // Empty sequences permitted
  DNA5SequenceCoding(DNA5SequenceCoding& copy) = delete; // For Performance reasons, no copy constructor.
  ~DNA5SequenceCoding() override = default;

  // For Performance reasons, don't allow naive assignments.
  DNA5SequenceCoding& operator=(DNA5SequenceCoding& copy) = delete;
  // Only allow move assignments
  DNA5SequenceCoding& operator=(DNA5SequenceCoding&& moved) noexcept {

    strand_ = moved.strand_;
    AlphabetSequence<CodingDNA5>::operator=(std::move(moved));
    return *this;

  }


  // Returns the sequence strand, FORWARD '+' or REVERSE '-'.
  [[nodiscard]] StrandSense strand() const { return strand_; }

private:

  StrandSense strand_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that cannot be used to directly generate an Amino Acid sequence
// This sequence is NOT STRANDED.
// However, it can return the STRANDED sequence DNA5SequenceCoding using a TranscriptionSequence (array of CDS).
// It can also be down-converted from a stranded sequence using the static linearSequence function.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A string of the standard 5 nucleotide DNA/RNA alphabet A, C, G, T/U, N
using StringDNA5 = AlphabetString<DNA5>;

// An UNSTRANDED DNA string.
class DNA5SequenceLinear: public  AlphabetSequence<DNA5> {

public:

  DNA5SequenceLinear(DNA5SequenceLinear&& sequence) noexcept : AlphabetSequence<DNA5>(std::move(sequence)) {}
  explicit DNA5SequenceLinear(StringDNA5&& sequence_string) :  AlphabetSequence<DNA5>(std::move(sequence_string)) {}
  explicit DNA5SequenceLinear(AlphabetSequence<DNA5>&& sequence) :  AlphabetSequence<DNA5>(std::move(sequence)) {}
  DNA5SequenceLinear() = default; // Empty sequences permitted.
  DNA5SequenceLinear(const DNA5SequenceLinear&) = delete; // For Performance reasons, no copy constructor.
  ~DNA5SequenceLinear() override = default;

  // For Performance reasons, don't allow naive assignments.
  DNA5SequenceLinear& operator=(const DNA5SequenceLinear&) = delete;
  // Only allow move assignments
  DNA5SequenceLinear& operator=(DNA5SequenceLinear&& moved) noexcept {

    AlphabetSequence<DNA5>::operator=(std::move(moved));
    return *this;

  }

  // Equality operator.
  [[nodiscard]] bool operator==(const DNA5SequenceLinear& cmp_seq) const { return equal(cmp_seq); }
  // Returns a defined subsequence (generally a single/group of codons) of the coding sequence
  // Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the TranscriptionFeatureMap.
  [[nodiscard]] DNA5SequenceCoding codingOffsetSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                            ContigOffset_t sub_sequence_offset,
                                                            ContigSize_t sub_sequence_length,
                                                            ContigOffset_t contig_offset) const;

  // Convenience routine that returns an array of sequences (strand adjusted), these will, in general, not be on codon boundaries.
  // Returned sequences are in transcription (strand) order with array[0] being the first exon.
  // The optional second offset argument is onlu used if the linear sequence is not a complete contig/chromosome.
  [[nodiscard]] std::vector<DNA5SequenceCoding> exonArraySequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                   ContigOffset_t contig_offset = 0) const;

  // The contig_offset adjusts for the offset in the contig from which the DNASequenceLinear was copied.
  // Setting sub_sequence_offset and sub_sequence_length to zero copies the entire intron sequence defined by the TranscriptionSequence.
  [[nodiscard]] DNA5SequenceCoding intronOffsetSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                            ContigOffset_t sub_sequence_offset,
                                                            ContigSize_t sub_sequence_length,
                                                            ContigOffset_t contig_offset) const;

  // Convenience routine that returns an array of introns (strand adjusted).
  // Returned sequences are in transcription (strand) order with array[0] being the first intron.
  // The optional second offset argument is only used if the linear sequence is not a complete contig/chromosome.
  [[nodiscard]] std::vector<DNA5SequenceCoding> intronArraySequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                     ContigOffset_t contig_offset = 0) const;

  // Offset is the relative sequence offset.
  [[nodiscard]] bool modifyBase(ContigOffset_t base_offset, DNA5::Alphabet Nucleotide);
  // Delete offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool deleteSubSequence(ContigOffset_t delete_offset, ContigSize_t delete_size);
  // Insert offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool insertSubSequence(ContigOffset_t insert_offset, const DNA5SequenceLinear& inserted_sequence);

  // Down-converts a coding DNA sequence to a linear DNA5 alphabet (swaps the logical alphabet from CodingDNA5 to DNA5).
  // Important - No strand conversion is performed.
  [[nodiscard]] static DNA5SequenceLinear downConvertToLinear(const DNA5SequenceCoding& stranded_sequence);


  // Down-converts a coding DNA sequence to a linear DNA5 alphabet (swaps the logical alphabet from CodingDNA5 to DNA5).
  // Important - Strand Conversion of reverse complement is performed for a -ve strand
  [[nodiscard]] static DNA5SequenceLinear strandedDownConvert(const DNA5SequenceCoding& stranded_sequence);

  // Up-converts a linear UNSTANDED DNA sequence to a STRANDED coding sequence (swaps the logical alphabet from DNA5 to CodingDNA5).
  // A -ve strand returns the reverse complement as expected.
  [[nodiscard]] DNA5SequenceCoding codingSequence(StrandSense strand) const;

  // Returns an UNSTRANDED subsequence. Returned sequence is valid but zero-sized if offset/size are out-of-bounds.
  [[nodiscard]] DNA5SequenceLinear subSequence(ContigOffset_t sub_sequence_offset, ContigSize_t sub_sequence_length) const;


private:

  // Down-converts a coding DNA sequence to a linear DNA5 alphabet (swaps the logical alphabet from CodingDNA5 to DNA5).
  // Important - Always returns the reverse complement of the coding sequence.
  [[nodiscard]] static DNA5SequenceLinear reverseDownConvert(const DNA5SequenceCoding& stranded_sequence);


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that implements a contig or chromosome.
// This sequence is NOT STRANDED.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UNSTRANDED
class DNA5SequenceContig: public DNA5SequenceLinear {

public:

  explicit DNA5SequenceContig(StringDNA5&& sequence_string) noexcept : DNA5SequenceLinear(std::move(sequence_string)) {}
  DNA5SequenceContig(const DNA5SequenceContig&) = delete; // For Performance reasons, no copy constructor.
  DNA5SequenceContig() = default;
  ~DNA5SequenceContig() override = default;

  // For Performance reasons, don't allow naive assignments.
  DNA5SequenceContig& operator=(const DNA5SequenceContig&) = delete;
  // Only allow move assignments
  DNA5SequenceContig& operator=(DNA5SequenceContig&& moved) noexcept {

    DNA5SequenceLinear::operator=(std::move(moved));
    return *this;

  }

  // Returns the codon offset of the contig offset within a coding sequence, returns false if not within the coding sequence.
  [[nodiscard]] bool codonOffset( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                  ContigOffset_t contig_offset,
                                  ContigOffset_t& codon_offset,
                                  ContigSize_t& base_in_codon) const;


  // Generally returns coding fragments such as Codons.
  [[nodiscard]] DNA5SequenceCoding codingSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                      ContigOffset_t sub_sequence_offset,
                                                      ContigSize_t sub_sequence_length) const {

    return codingOffsetSubSequence(coding_seq_ptr, sub_sequence_offset, sub_sequence_length, 0);

  }

  // The entire sequence defined by the TranscriptionSequence is returned.
  [[nodiscard]] DNA5SequenceCoding codingSequence(const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr) const {

    return codingSubSequence(coding_seq_ptr, 0, 0);

  }

  // The entire intron sequence defined by the TranscriptionSequence is returned.
  [[nodiscard]] DNA5SequenceCoding intronSequence(const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr) const {

    return intronOffsetSubSequence(coding_seq_ptr, 0, 0, 0);

  }


private:


};



}   // end namespace


#endif //KGL_SEQUENCE_BASE_H
