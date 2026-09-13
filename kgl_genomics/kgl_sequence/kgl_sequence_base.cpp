//
// Created by kellerberrin on 31/10/17.
//


#include "kgl_sequence_base.h"
#include "kel_interval_unsigned.h"

#include <ranges>

namespace kgl = kellerberrin::genome;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that cannot be used to directly generate an Amino Acid sequence
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The base DNA5 sequence class.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Letter offset is relative to the begining of the sequence (0 is the first letter).
bool kgl::DNA5SequenceLinear::modifyBase(ContigOffset_t base_offset, DNA5::Alphabet nucleotide) {

  return modifyLetter(base_offset, nucleotide);

}

bool kgl::DNA5SequenceLinear::deleteSubSequence(const OpenRightUnsigned& delete_interval) {

  return deleteOffset(delete_interval.lower(), delete_interval.size());

}


// Insert offset is relative to the begining of the sequence (0 is the first letter).
bool kgl::DNA5SequenceLinear::insertSubSequence(ContigOffset_t insert_offset, const DNA5SequenceLinear& inserted_sequence) {

  return insertOffset(insert_offset, inserted_sequence);

}


// Returns an UNSTRANDED subsequence. Returned sequence is valid but zero-sized if offset/size are out-of-bounds.
std::optional<kgl::DNA5SequenceLinear> kgl::DNA5SequenceLinear::subSequence(const OpenRightUnsigned& sub_interval) const {

  auto sub_sequence_opt = getSubsequence(sub_interval);
  if (not sub_sequence_opt) {

    ExecEnv::log().error("Cannot get sub-sequence: {} from sequence: {}", sub_interval.toString(), interval().toString());
    return std::nullopt;

  }

  return DNA5SequenceLinear(std::move(sub_sequence_opt.value()));

}


// Down-converts a coding DNA sequence to a linear DNA5 alphabet (swaps the logical alphabet from CodingDNA5 to DNA5).
// Important - No strand conversion is performed.
kgl::DNA5SequenceLinear kgl::DNA5SequenceLinear::downConvertToLinear(const DNA5SequenceCoding& stranded_sequence) {

  StringDNA5 linear_string;
  linear_string.reserve(stranded_sequence.length()); // pre-allocate for efficiency

  auto convert_base = [](CodingDNA5::Alphabet base)->DNA5::Alphabet { return DNA5::convertFromCodingDNA5(base); };
  std::ranges::transform(stranded_sequence.getAlphabetString(), std::back_inserter(linear_string), convert_base);

  return DNA5SequenceLinear(std::move(linear_string));

}


// Up-converts a linear UNSTANDED DNA sequence to a STRANDED coding sequence (swaps the logical alphabet from DNA5 to CodingDNA5).
// A -ve strand returns the reverse complement as expected.
// Delegates to the sequence view, which operates on exactly the same bytes with the same transforms.
kgl::DNA5SequenceCoding kgl::DNA5SequenceLinear::codingSequence(StrandSense strand) const {

  return getView().codingSequence(strand);

}

std::optional<kgl::DNA5SequenceLinear> kgl::DNA5SequenceLinear::concatSequences(const IntervalSetLower& interval_set) const {

  // Concatenate the extracted sub-sequence views in interval (sorted) order.
  // The first sub-view initialises the result; the rest are appended.
  std::optional<DNA5SequenceLinear> concatenated_sequence;
  for (auto const& sub_interval : interval_set) {

    auto sub_view_opt = getView().subView(sub_interval);
    if (not sub_view_opt) {

      ExecEnv::log().warn("Unable to extract sub-sequence: {} for interval: {}", sub_interval.toString(), interval().toString());
      return std::nullopt;

    }

    if (concatenated_sequence) {

      if (not concatenated_sequence->append(DNA5SequenceLinear(sub_view_opt.value()))) {

        ExecEnv::log().warn("Unable to concatenate modified sequence for interval");
        return std::nullopt;

      }

    } else {

      concatenated_sequence.emplace(sub_view_opt.value());

    }

  }

  if (not concatenated_sequence) {

    ExecEnv::log().warn("No concat sub-sequences for interval set size: {}", interval_set.size());
    return std::nullopt;

  }

  return concatenated_sequence;

}