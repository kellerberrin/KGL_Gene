//
// Created by kellerberrin on 20/11/23.
//

#ifndef KGL_SEQUENCE_BASE_VIEW_H
#define KGL_SEQUENCE_BASE_VIEW_H

#include "kgl_alphabet_view.h"
#include "kgl_genome_prelim.h"

#include <string>
#include <memory>

namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Amino Sequence View - A view on Amino Acid (protein) sequences.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward.
class AminoSequence;

class AminoSequenceView: public AlphabetView<AminoAcid> {

public:

  explicit AminoSequenceView(const AminoSequence& sequence);
  AminoSequenceView(const AminoSequenceView& amino_view) = default;
  ~AminoSequenceView() override = default;

  [[nodiscard]] std::optional<AminoSequenceView> subView(const OpenRightUnsigned& sub_interval) const;

private:

  explicit AminoSequenceView(const AlphabetView<AminoAcid>& alphabet_view) : AlphabetView<AminoAcid>(alphabet_view) {}

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A view on a STRANDED DNA string that can be converted to an AMINO sequence.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward.
class DNA5SequenceCoding;

class DNA5SequenceCodingView: public AlphabetView<CodingDNA5> {

public:

  explicit DNA5SequenceCodingView(const DNA5SequenceCoding& coding_sequence);
  DNA5SequenceCodingView(const DNA5SequenceCodingView& copy) = default;
  ~DNA5SequenceCodingView() override = default;

  [[nodiscard]] std::optional<DNA5SequenceCodingView> subView(const OpenRightUnsigned& sub_interval) const;

  // Returns the sequence strand, FORWARD '+' or REVERSE '-'.
  [[nodiscard]] StrandSense strand() const { return strand_; }

private:

  StrandSense strand_;

  DNA5SequenceCodingView(const AlphabetView<CodingDNA5>& sequence, StrandSense strand) : AlphabetView<CodingDNA5>(sequence), strand_(strand) {}

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A view of a  linear and contiguous DNA5 sequence that cannot be used to directly generate an Amino Acid sequence
// This sequence is NOT STRANDED.
// However, it can return the STRANDED sequence DNA5SequenceCoding using a TranscriptionSequence (array of CDS).
// It can also be down-converted from a stranded sequence using the static linearSequence function.
// A string of the standard 5 nucleotide DNA/RNA alphabet A, C, G, T/U, N
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward.
class DNA5SequenceLinear;

class DNA5SequenceLinearView: public  AlphabetView<DNA5> {

public:

  explicit DNA5SequenceLinearView(const DNA5SequenceLinear& sequence);
  DNA5SequenceLinearView(const DNA5SequenceLinearView& copy) = default;
  ~DNA5SequenceLinearView() override = default;

  // Up-converts a linear UNSTANDED DNA sequence to a STRANDED coding sequence (swaps the logical alphabet from DNA5 to CodingDNA5).
  // A -ve strand returns the reverse complement as expected.
  [[nodiscard]] DNA5SequenceCoding codingSequence(StrandSense strand) const;
  [[nodiscard]] std::optional<DNA5SequenceLinearView> subView(const OpenRightUnsigned& sub_interval) const;
  
private:

  explicit DNA5SequenceLinearView(const AlphabetView<DNA5>& sequence) :  AlphabetView<DNA5>(sequence) {}

};



}   // end namespace

#endif //KGL_SEQUENCE_BASE_VIEW_H
