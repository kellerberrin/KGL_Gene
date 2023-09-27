//
// Created by kellerberrin on 31/10/17.
//


#include "kgl_sequence_base.h"
#include "kgl_genome_seq/kgl_seq_interval.h"
#include "kel_interval_set.h"

#include <ranges>

namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that cannot be used to directly generate an Amino Acid sequence
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The base DNA5 sequence class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
std::optional<kgl::DNA5SequenceLinear> kgl::DNA5SequenceLinear::subOptSequence(const OpenRightUnsigned& sub_interval) const {

  if (not interval().containsInterval(sub_interval)) {

    ExecEnv::log().warn("Sub interval: {} not contained in interval: {}.", sub_interval.toString(), interval().toString());
    return std::nullopt;

  }

  DNA5SequenceLinear sub_sequence;

  if (not getSubsequence(sub_interval.lower(), sub_interval.size(), sub_sequence)) {

    ExecEnv::log().warn("Cannot get sub-sequence: {} from sequence: {}", sub_interval.toString(), interval().toString());
    // ReturnType an empty sequence
    return sub_sequence;

  }

  return sub_sequence;

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
kgl::DNA5SequenceCoding kgl::DNA5SequenceLinear::codingSequence(StrandSense strand) const {


  StringCodingDNA5 coding_string;
  coding_string.reserve(length()); // pre-allocate for efficiency.

  if (strand == StrandSense::REVERSE) {

    auto complement_coding = [](DNA5::Alphabet base)->CodingDNA5::Alphabet { return DNA5::complementNucleotide(base); };
    auto reversed_string = std::ranges::views::reverse(getAlphabetString());
    std::ranges::transform(reversed_string, std::back_inserter(coding_string), complement_coding);

  } else { // Strand is positive or unknown.

    auto convert_coding = [](DNA5::Alphabet base)->CodingDNA5::Alphabet { return DNA5::convertToCodingDNA5(base); };
    std::ranges::transform(getAlphabetString(), std::back_inserter(coding_string), convert_coding);

  }

  return DNA5SequenceCoding(std::move(coding_string), strand);

}

std::optional<kgl::DNA5SequenceLinear> kgl::DNA5SequenceLinear::concatSequences(const std::vector<OpenRightUnsigned>& interval_vector) const {

  // Sort the intervals.
  IntervalSetLower interval_set;
  for (auto const& interval : interval_vector) {

    auto [insert_iter, result] = interval_set.insert(interval);
    if (not result) {

      ExecEnv::log().warn("Interval: {} has a duplicate lower()", interval.toString());

    }

  }

  // Extract the modified sequences and concatenate them.
  DNA5SequenceLinear concatenated_sequence;
  for (auto const& sub_interval : interval_set) {

    if (not interval().containsInterval(sub_interval)) {

      ExecEnv::log().warn("Sub-interval: {} unable is not contained in interval: {}",
                          sub_interval.toString(), interval().toString());
      return std::nullopt;

    }

    auto sub_sequence_opt = subOptSequence(sub_interval);
    if (not sub_sequence_opt) {

      ExecEnv::log().warn("Unable to extract sub-sequence: {} for interval: {}", sub_interval.toString(), interval().toString());
      return sub_sequence_opt;

    }

    bool result = concatenated_sequence.append(sub_sequence_opt.value());
    if (not result) {

      ExecEnv::log().warn("Unable to concatenate modified sequence for interval: {}",
                          sub_interval.toString());
      return std::nullopt;

    }

  }

  return concatenated_sequence;

}

