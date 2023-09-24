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

// Delete offset is relative to the begining of the sequence (0 is the first letter).
bool kgl::DNA5SequenceLinear::deleteSubSequence(ContigOffset_t delete_offset, ContigSize_t delete_size) {

  return deleteOffset(delete_offset, delete_size);

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

    ExecEnv::log().warn_location("Sub interval: {} not contained in interval: {}.", sub_interval.toString(), interval().toString());
    return std::nullopt;

  }

  return subSequence(sub_interval.lower(), sub_interval.size());

}



// Returns an UNSTRANDED subsequence. Returned sequence is valid but zero-sized if offset/size are out-of-bounds.
kgl::DNA5SequenceLinear kgl::DNA5SequenceLinear::subSequence(const OpenRightUnsigned& sub_interval) const {

  return subSequence(sub_interval.lower(), sub_interval.size());

}

kgl::DNA5SequenceLinear kgl::DNA5SequenceLinear::subSequence(ContigOffset_t offset, ContigSize_t sub_length) const {

  DNA5SequenceLinear sub_sequence;

  // Check offset and size.
  if ((offset + sub_length) > length() or sub_length > length()) {

    ExecEnv::log().warn("DNA5SequenceLinear::subSequence; sub sequence offset: {} and sub sequence size: {} too large for sequence length: {}",
                        offset, sub_length, length());
    // ReturnType an empty sequence
    return sub_sequence;

  }

  if (not getSubsequence(offset, sub_length, sub_sequence)) {

    ExecEnv::log().warn("DNA5SequenceLinear::subSequence; Cannot get sub-sequence offset: {} and sub sequence size: {} from sequence length: {}",
                        offset, sub_length, length());
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

// Down-converts a coding DNA sequence to a linear DNA5 alphabet (swaps the logical alphabet from CodingDNA5 to DNA5).
// Important - Always returns the reverse complement of the stranded sequence.
kgl::DNA5SequenceLinear kgl::DNA5SequenceLinear::reverseDownConvert(const DNA5SequenceCoding& stranded_sequence) {

  StringDNA5 linear_string;
  linear_string.reserve(stranded_sequence.length()); // pre-allocate for efficiency

  auto complement = [](CodingDNA5::Alphabet coding_base)->DNA5::Alphabet { return DNA5::convertComplementNucleotide(coding_base); };
  auto reversed_string = std::ranges::views::reverse(stranded_sequence.getAlphabetString());
  std::ranges::transform(reversed_string,std::back_inserter(linear_string), complement);

  return DNA5SequenceLinear(std::move(linear_string));

}

// Down-converts a coding DNA sequence to a linear DNA5 alphabet (swaps the logical alphabet from CodingDNA5 to DNA5).
// Important - Strand Conversion of a reverse complement is performed for a -ve strand
kgl::DNA5SequenceLinear kgl::DNA5SequenceLinear::strandedDownConvert(const DNA5SequenceCoding& stranded_sequence) {

  switch(stranded_sequence.strand()) {

    case StrandSense::FORWARD:
      return downConvertToLinear(stranded_sequence);

    case StrandSense::REVERSE:
      return reverseDownConvert(stranded_sequence);

  }

  return downConvertToLinear(stranded_sequence); // To keep the compiler happy

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

      ExecEnv::log().warn("DNA5SequenceLinear::concatSequences; interval: {} has a duplicate lower()", interval.toString());

    }

  }

  // Extract the modified sequences and concatenate them.
  DNA5SequenceLinear concatenated_sequence;
  for (auto const& sub_interval : interval_set) {

    if (not interval().containsInterval(sub_interval)) {

      ExecEnv::log().warn("DNA5SequenceLinear::concatSequences; sub-interval: {} unable is not contained in interval: {}",
                          sub_interval.toString(), interval().toString());
      return std::nullopt;

    }

    auto sub_sequence = subSequence(sub_interval);

    bool result = concatenated_sequence.append(sub_sequence);
    if (not result) {

      ExecEnv::log().warn("SequenceTranscript::concatOriginalSequences; unable to concatenate modified sequence for interval: {}",
                          sub_interval.toString());
      return std::nullopt;

    }

  }

  return concatenated_sequence;

}


kgl::DNA5SequenceCoding kgl::DNA5SequenceLinear::codingSequence(const std::shared_ptr<const TranscriptionSequence>& transcript_ptr) const {

  auto cds_interval_set = GeneIntervalStructure::transcriptIntervals(transcript_ptr);
  std::vector<OpenRightUnsigned> interval_vector(cds_interval_set.begin(), cds_interval_set.end());

  auto concat_sequence_opt = concatSequences(interval_vector);
  if (not concat_sequence_opt) {

    ExecEnv::log().warn("DNA5SequenceLinear::codingSequence; Unable to concat sequence intervals for Gene: {}, Transcript: {}",
                        transcript_ptr->getGene()->id(), transcript_ptr->getParent()->id());

    return {}; // Return an empty coding sequence.

  }

  auto& concat_sequence = concat_sequence_opt.value();

  return concat_sequence.codingSequence(transcript_ptr->strand());

}

