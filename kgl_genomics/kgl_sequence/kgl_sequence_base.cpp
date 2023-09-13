//
// Created by kellerberrin on 31/10/17.
//


#include "kgl_sequence_base.h"
#include "kgl_sequence_codon.h"
#include "kgl_sequence_offset.h"

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


// The contig_offset adjusts for the offset in the contig from which the DNASequenceLinear was copied.
// Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the TranscriptionSequence.
kgl::DNA5SequenceCoding kgl::DNA5SequenceLinear::codingOffsetSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                          ContigOffset_t sub_sequence_offset,
                                                                          ContigSize_t sub_sequence_length,
                                                                          ContigOffset_t contig_offset) const {

  return SequenceOffset::refCodingSubSequence(coding_seq_ptr, *this, sub_sequence_offset, sub_sequence_length, contig_offset);

}


// The contig_offset adjusts for the offset in the contig from which the DNASequenceLinear was copied.
// Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the TranscriptionSequence.
std::vector<kgl::DNA5SequenceCoding> kgl::DNA5SequenceLinear::exonArraySequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                                 ContigOffset_t contig_offset) const {

  return SequenceOffset::refExonArraySequence(coding_seq_ptr, *this, 0, 0, contig_offset);

}


// The contig_offset adjusts for the offset in the contig from which the DNASequenceLinear was copied.
// Setting sub_sequence_offset and sub_sequence_length to zero copies the entire intron sequence defined by the TranscriptionSequence.
kgl::DNA5SequenceCoding kgl::DNA5SequenceLinear::intronOffsetSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                          ContigOffset_t sub_sequence_offset,
                                                                          ContigSize_t sub_sequence_length,
                                                                          ContigOffset_t contig_offset) const {

  return SequenceOffset::refIntronSubSequence(coding_seq_ptr, *this, sub_sequence_offset, sub_sequence_length, contig_offset);

}

// Convenience routine that returns an array of introns (strand adjusted).
// Returned sequences are in transcription (strand) order with array[0] being the first intron.
// The optional second offset argument is onlu used if the linear sequence is not a complete contig/chromosome.
std::vector<kgl::DNA5SequenceCoding> kgl::DNA5SequenceLinear::intronArraySequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                                   ContigOffset_t contig_offset) const {

  return SequenceOffset::refIntronArraySequence(coding_seq_ptr, *this, 0, 0, contig_offset);

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

  auto convert_base = [](CodingDNA5::Alphabet base) { return DNA5::convertFromCodingDNA5(base); };
  std::transform(stranded_sequence.getAlphabetString().begin(),
                 stranded_sequence.getAlphabetString().end(),
                 std::back_inserter(linear_string), convert_base);

  return DNA5SequenceLinear(std::move(linear_string));

}

// Down-converts a coding DNA sequence to a linear DNA5 alphabet (swaps the logical alphabet from CodingDNA5 to DNA5).
// Important - Always returns the reverse complement of the stranded sequence.
kgl::DNA5SequenceLinear kgl::DNA5SequenceLinear::reverseDownConvert(const DNA5SequenceCoding& stranded_sequence) {

  StringDNA5 linear_string;
  linear_string.reserve(stranded_sequence.length()); // pre-allocate for efficiency

  auto complement_base = [](CodingDNA5::Alphabet coding_base) { return DNA5::convertComplementNucleotide(coding_base); };
  std::transform(stranded_sequence.getAlphabetString().rbegin(),
                 stranded_sequence.getAlphabetString().rend(),
                 std::back_inserter(linear_string), complement_base);

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
  coding_string.reserve(length()); // pre-allocate for efficiency

  if (strand == StrandSense::REVERSE) {

    for (auto rit = getAlphabetString().rbegin(); rit != getAlphabetString().rend(); ++rit) {

      coding_string.push_back(DNA5::complementNucleotide(*rit));

    }

  } else { // Strand is positive or unknown.

    auto convert_coding = [](DNA5::Alphabet base) { return DNA5::convertToCodingDNA5(base); };

    std::transform(getAlphabetString().begin(),
                   getAlphabetString().end(),
                   std::back_inserter(coding_string),
                   convert_coding);

  }

  return DNA5SequenceCoding(std::move(coding_string), strand);

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A STRANDED DNA string that can be converted to an AMINO sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence used in a contig (chromosome). This object exists for semantic reasons.
// NOT STRANDED
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Returns the codon offset of offset within a coding, returns false if not within the coding sequence.
bool kgl::DNA5SequenceContig::codonOffset(const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                          ContigOffset_t contig_offset,
                                          ContigOffset_t& codon_offset,
                                          ContigSize_t& base_in_codon) const {

  ContigOffset_t sequence_offset;
  ContigSize_t sequence_length;
  if (SequenceOffset::refOffsetWithinCodingSequence(coding_seq_ptr, contig_offset, sequence_offset, sequence_length)) {

    codon_offset = static_cast<ContigOffset_t>(sequence_offset / Codon::CODON_SIZE);
    base_in_codon = static_cast <ContigOffset_t>(sequence_offset % Codon::CODON_SIZE);
    return true;

  } else {

    codon_offset = 0;
    base_in_codon = 0;
    return false;

  }

}

