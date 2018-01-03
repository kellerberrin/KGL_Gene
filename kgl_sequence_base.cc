//
// Created by kellerberrin on 31/10/17.
//


#include "kgl_sequence_base.h"
#include "kgl_sequence_codon.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The base DNA5 sequence class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Letter offset is relative to the begining of the sequence (0 is the first letter).
bool kgl::DNA5SequenceCoding::modifyBase(ContigOffset_t sequence_offset, CodingDNA5::Alphabet nucleotide) {

  return modifyLetter(sequence_offset, nucleotide);

}

// Delete offset is relative to the begining of the sequence (0 is the first letter).
bool kgl::DNA5SequenceCoding::deleteSubSequence(ContigOffset_t delete_offset, ContigSize_t delete_size) {

  return deleteOffset(delete_offset, delete_size);

}

// Insert offset is relative to the begining of the sequence (0 is the first letter).
bool kgl::DNA5SequenceCoding::insertSubSequence(ContigOffset_t insert_offset, const DNA5SequenceCoding& inserted_sequence) {

  return insertOffset(insert_offset, inserted_sequence);

}




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

// Insert offset is relative to the begining of the sequence (0 is the first letter).
bool kgl::DNA5SequenceLinear::insertSubSequence(ContigOffset_t insert_offset, const DNA5SequenceLinear& inserted_sequence) {

  return insertOffset(insert_offset, inserted_sequence);

}


// The contig_offset adjusts for the offset in the contig from which the DNASequenceLinear was copied.
// Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the CodingSequence.
std::shared_ptr<kgl::DNA5SequenceCoding>
kgl::DNA5SequenceLinear::codingOffsetSubSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                 ContigOffset_t sub_sequence_offset,
                                                 ContigSize_t sub_sequence_length,
                                                 ContigOffset_t contig_offset) const {

  return SequenceOffset::refCodingSubSequence(coding_seq_ptr, *this, sub_sequence_offset, sub_sequence_length, contig_offset);

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence used in a contig (chromosome). This object exists for semantic reasons.
// NOT STRANDED
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Returns the codon offset of offset within a coding, returns false if not within the coding sequence.
bool kgl::DNA5SequenceContig::codonOffset(std::shared_ptr<const CodingSequence> coding_seq_ptr,
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

