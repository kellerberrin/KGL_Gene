//
// Created by kellerberrin on 3/01/18.
//

#ifndef KGL_SEQUENCE_OFFSET_H
#define KGL_SEQUENCE_OFFSET_H

#include "kgl_sequence_base.h"
#include "kgl_variant_mutation_offset.h"



namespace kellerberrin::genome {   //  organization level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A sorted map of coding sequence offsets.
// The first (sort key) offset is the beginning of the EXON and the second (value) is the EXON end offset.
// The offsets use the half interval idiom [start, end).
// The offsets are initially contig offsets. However these can be modified to account for indel modifications.
// The map is generated by an adapter function exonOffsetAdapter() that uses a coding sequence as input.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using ExonOffsetMap = std::map<ContigOffset_t, ContigOffset_t>;
using IntronOffsetMap = ExonOffsetMap;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This class implements the tricky offset logic associated with generating a coding sequence after the underlying
// DNA sequence has been modified by indel variants.
// A Singleton class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class SequenceOffset {

public:

  SequenceOffset() = delete;
  ~SequenceOffset() = delete;

// Convenience routine that returns a coding sequence from an unmutated (reference) sequence
  [[nodiscard]] static DNA5SequenceCoding refCodingSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                const DNA5SequenceLinear& sequence,
                                                                ContigOffset_t sub_sequence_offset,
                                                                ContigSize_t sub_sequence_length,
                                                                ContigOffset_t contig_offset);

// Convenience routine that returns an array of sequences (strand adjusted) from an unmutated (reference) sequence.
// Returned sequences are in transcription (strand) order with array[0] being the first exon.
  [[nodiscard]] static std::vector<DNA5SequenceCoding> refExonArraySequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                             const DNA5SequenceLinear& sequence,
                                                                             ContigOffset_t sub_sequence_offset,
                                                                             ContigSize_t sub_sequence_length,
                                                                             ContigOffset_t contig_offset);


// Returns a defined subsequence of all the introns of the coding sequence concatonated.
// Setting sub_sequence_offset and sub_sequence_length to zero copies the entire intron sequence defined by the TranscriptionFeatureMap.
  [[nodiscard]] static DNA5SequenceCoding refIntronSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                const DNA5SequenceLinear& sequence_ptr,
                                                                ContigOffset_t sub_sequence_offset,
                                                                ContigSize_t sub_sequence_length,
                                                                ContigOffset_t contig_offset);

// Convenience routine that returns an array of introns (strand adjusted) from an unmutated (reference) sequence
// Returned sequences are in transcription (strand) order with array[0] being the first intron.
  [[nodiscard]] static std::vector<DNA5SequenceCoding> refIntronArraySequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                               const DNA5SequenceLinear& sequence,
                                                                               ContigOffset_t sub_sequence_offset,
                                                                               ContigSize_t sub_sequence_length,
                                                                               ContigOffset_t contig_offset);

  // Convenience routine that returns a coding sequence from an unmutated (reference) sequence
  [[nodiscard]] static DNA5SequenceCoding mutantCodingSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                   const DNA5SequenceLinear& sequence,
                                                                   const VariantMutationOffset& indel_adjust,
                                                                   ContigOffset_t sub_sequence_offset,
                                                                   ContigSize_t sub_sequence_length,
                                                                   ContigOffset_t contig_offset);

// Returns bool false if contig_offset is not within the coding sequence defined by the coding_seq_ptr.
// If the contig_offset is in the coding sequence then a valid sequence_offset and the sequence length is returned.
// The offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
  [[nodiscard]] static bool refOffsetWithinCodingSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                           ContigOffset_t contig_offset,
                                                           ContigOffset_t &coding_sequence_offset,
                                                           ContigSize_t &coding_sequence_length);


    // Inverse of the above function. Given the stranded base offset within a coding sequence, return the corresponding contig offset.
// Returns bool false if then coding sequence offset is not within the coding sequence defined by the coding_seq_ptr.
// If the coding sequence_offset is in the coding sequence then a valid contig_offset and the sequence length is returned.
// The contig offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
  [[nodiscard]] static bool refCodingSequenceContigOffset( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                           ContigOffset_t coding_sequence_offset,
                                                           ContigOffset_t &contig_offset,
                                                           ContigSize_t &coding_sequence_length);

// A sorted map of intron coding sequence offsets.
// The first (sort key) offset is the beginning of the intron and the second (value) is the intron end offset.
// The offsets use the half interval idiom [start, end).
// The offsets are initially contig offsets. However these can be modified to account for indel offset modification.
// The map is generated by an adapter function intronOffsetAdapter() that uses a coding sequence as input.
// If the coding sequence only defines 1 exon then an empty map is returned.
// In general the number of intron offsets returned is the number of defined sequences - 1
  [[nodiscard]] static bool intronOffsetAdapter( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                 StrandSense& strand,
                                                 IntronOffsetMap& intron_offset_map);

  // Converts linear DNA to a coding DNA sequence.
  [[nodiscard]] static DNA5SequenceCoding codingSequence(const DNA5SequenceLinear& base_sequence, StrandSense strand);

private:

  // Returns a defined subsequence (generally a single/group of codons) of the coding sequence
  // Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the TranscriptionFeatureMap.
  [[nodiscard]] static DNA5SequenceCoding codingSubSequence( const DNA5SequenceLinear& base_sequence,
                                                             const ExonOffsetMap& exon_offset_map,
                                                             StrandSense strand,
                                                             ContigOffset_t sub_sequence_offset, // base count
                                                             ContigSize_t sub_sequence_length,  // number of bases
                                                             ContigOffset_t contig_offset); // if a subsequence.

// A sorted map of coding sequence offsets.
// The first (sort key) offset is the beginning of the EXON and the second (value) is the EXON end offset.
// The offsets use the half interval idiom [start, end).
// The offsets are initially contig offsets. However these can be modified to account for indel offset modification.
// The map is generated by an adapter function exonOffsetAdapter() that uses a coding sequence as input.
  [[nodiscard]] static bool exonOffsetAdapter( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                               StrandSense& strand,
                                               ExonOffsetMap& exon_offset_map);

  // Adjusted for indel mutations.
  [[nodiscard]] static bool exonMutantOffset( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                              const VariantMutationOffset& indel_adjust,
                                              StrandSense& strand,
                                              ExonOffsetMap& exon_offset_map);


// Returns bool false if contig_offset is not within the coding sequence defined by the coding_seq_ptr.
// If the contig_offset is in the coding sequence then a valid sequence_offset and the sequence length is returned.
// The offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
  [[nodiscard]] static bool offsetWithinCodingSequence( const ExonOffsetMap& exon_offset_map,
                                                        StrandSense strand,
                                                        ContigOffset_t sequence_offset,
                                                        ContigOffset_t sequence_contig_offset,
                                                        ContigOffset_t &coding_sequence_offset,
                                                        ContigSize_t &coding_sequence_length);

// Inverse of the above function. Returns bool false if the sequence_offset is not within the coding sequence defined by the coding_seq_ptr.
// If the sequence_offset is in the coding sequence then a valid contig_offset and the sequence length is returned.
// The contig offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
  [[nodiscard]] static bool codingSequenceContigOffset( const ExonOffsetMap& exon_offset_map,
                                                        StrandSense strand,
                                                        ContigOffset_t sequence_offset,
                                                        ContigOffset_t &contig_offset,
                                                        ContigSize_t &coding_sequence_length);

};



}   // end namespace


#endif //KGL_SEQUENCE_OFFSET_H
