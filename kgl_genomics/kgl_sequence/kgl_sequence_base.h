//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SEQUENCE_BASE_H
#define KGL_SEQUENCE_BASE_H


#include "kgl_alphabet_string.h"
#include "kgl_alphabet_sequence.h"
#include "kgl_genome_prelim.h"
#include "kgl_sequence_base_view.h"

#include <string>
#include <memory>

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
// This sequence is ALWAYS in strand sense (reverse complement has already been performed on -ve strand sequences).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A string of the standard 5 nucleotide DNA/RNA alphabet A, C, G, T/U, N
using StringCodingDNA5 = AlphabetString<CodingDNA5>;

// A STRANDED DNA string that can be converted to an AMINO sequence.
class DNA5SequenceCoding: public AlphabetSequence<CodingDNA5> {

public:

  DNA5SequenceCoding(DNA5SequenceCoding&& sequence) noexcept :  AlphabetSequence<CodingDNA5>(std::move(sequence.alphabet_string_)), strand_(sequence.strand_) {}
  explicit DNA5SequenceCoding(const DNA5SequenceCodingView& sequence_view) : AlphabetSequence<CodingDNA5>(sequence_view.getSequence()), strand_(sequence_view.strand()) {};
  DNA5SequenceCoding(StringCodingDNA5&& sequence_string, StrandSense strand) noexcept : AlphabetSequence<CodingDNA5>(std::move(sequence_string)), strand_(strand) {}
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

  //Return a view.
  [[nodiscard]] DNA5SequenceCodingView getView() const { return DNA5SequenceCodingView(*this); }

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
  explicit DNA5SequenceLinear(const DNA5SequenceLinearView& sequence_view) : AlphabetSequence<DNA5>(sequence_view.getSequence()) {};
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

  // Offset is the relative sequence offset.
  [[nodiscard]] bool modifyBase(ContigOffset_t base_offset, DNA5::Alphabet Nucleotide);
  // Delete offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool deleteSubSequence(const OpenRightUnsigned& delete_interval);
  // Insert offset is relative to the begining of the sequence (0 is the first letter).
  [[nodiscard]] bool insertSubSequence(ContigOffset_t insert_offset, const DNA5SequenceLinear& inserted_sequence);

  // Down-converts a coding DNA sequence to a linear DNA5 alphabet (swaps the logical alphabet from CodingDNA5 to DNA5).
  // Important - No strand conversion is performed.
  [[nodiscard]] static DNA5SequenceLinear downConvertToLinear(const DNA5SequenceCoding& stranded_sequence);

  // Up-converts a linear UNSTANDED DNA sequence to a STRANDED coding sequence (swaps the logical alphabet from DNA5 to CodingDNA5).
  // A -ve strand returns the reverse complement as expected.
  [[nodiscard]] DNA5SequenceCoding codingSequence(StrandSense strand) const;

  // Returns an UNSTRANDED subsequence.
  [[nodiscard]] std::optional<DNA5SequenceLinear> subSequence(const OpenRightUnsigned& sub_interval) const;
  //Return a view.
  [[nodiscard]] DNA5SequenceLinearView getView() const { return DNA5SequenceLinearView(*this); }

  // Sorts a vector of intervals in lower() ascending order and then concatanates the sub-intervals together.
  // Does not check for overlapping intervals.
  [[nodiscard]] std::optional<DNA5SequenceLinear> concatSequences(const IntervalSetLower& interval_set) const;

private:


};



}   // end namespace


#endif //KGL_SEQUENCE_BASE_H
