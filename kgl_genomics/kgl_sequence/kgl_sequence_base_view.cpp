//
// Created by kellerberrin on 23/11/23.
//

#include "kgl_sequence_amino.h"
#include "kgl_sequence_base.h"
#include "kgl_sequence_base_view.h"


#include <ranges>


namespace kgl = kellerberrin::genome;


kgl::AminoSequenceView::AminoSequenceView(const AminoSequence& sequence)
: AlphabetView<AminoAcid>(sequence) {}


kgl::DNA5SequenceCodingView::DNA5SequenceCodingView(const DNA5SequenceCoding& coding_sequence)
: AlphabetView<CodingDNA5>(coding_sequence), strand_(coding_sequence.strand()) {}


kgl::DNA5SequenceLinearView::DNA5SequenceLinearView(const DNA5SequenceLinear& sequence)
: AlphabetView<DNA5>(sequence) {}


// Returns subsequence view. getSubView() validates the interval and logs a warning if out of bounds.
std::optional<kgl::AminoSequenceView> kgl::AminoSequenceView::subView(const OpenRightUnsigned& sub_interval) const {

  auto sub_view_opt = getSubView(sub_interval);
  if (not sub_view_opt) return std::nullopt;
  return AminoSequenceView(sub_view_opt.value());

}

// Returns subsequence view. getSubView() validates the interval and logs a warning if out of bounds.
std::optional<kgl::DNA5SequenceCodingView> kgl::DNA5SequenceCodingView::subView(const OpenRightUnsigned& sub_interval) const {

  auto sub_view_opt = getSubView(sub_interval);
  if (not sub_view_opt) return std::nullopt;
  return DNA5SequenceCodingView(sub_view_opt.value(), strand());

}


// Returns subsequence view. getSubView() validates the interval and logs a warning if out of bounds.
std::optional<kgl::DNA5SequenceLinearView> kgl::DNA5SequenceLinearView::subView(const OpenRightUnsigned& sub_interval) const {

  auto sub_view_opt = getSubView(sub_interval);
  if (not sub_view_opt) return std::nullopt;
  return DNA5SequenceLinearView(sub_view_opt.value());

}


// Up-converts a linear UNSTANDED DNA sequence to a STRANDED coding sequence (swaps the logical alphabet from DNA5 to CodingDNA5).
// A -ve strand returns the reverse complement as expected.
kgl::DNA5SequenceCoding kgl::DNA5SequenceLinearView::codingSequence(StrandSense strand) const {


  StringCodingDNA5 coding_string;
  coding_string.reserve(length()); // pre-allocate for efficiency.

  if (strand == StrandSense::REVERSE) {

    auto complement_coding = [](DNA5::Alphabet base)->CodingDNA5::Alphabet { return DNA5::complementNucleotide(base); };
    auto reversed_string = std::ranges::views::reverse(alphabet_view_);
    std::ranges::transform(reversed_string, std::back_inserter(coding_string), complement_coding);

  } else {

    auto convert_coding = [](DNA5::Alphabet base)->CodingDNA5::Alphabet { return DNA5::convertToCodingDNA5(base); };
    std::ranges::transform(alphabet_view_, std::back_inserter(coding_string), convert_coding);

  }

  return DNA5SequenceCoding(std::move(coding_string), strand);

}