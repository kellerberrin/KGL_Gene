//
// Created by kellerberrin on 23/09/23.
//

// Remove ASAP

#include "kgl_sequence_base.h"

namespace kgl = kellerberrin::genome;

kgl::DNA5SequenceCoding
kgl::DNA5SequenceLinear::intronSequence( const std::shared_ptr<const TranscriptionSequence>&,
                                         ContigOffset_t) const { return DNA5SequenceCoding(); }

// Convenience routine that returns an array of introns (strand adjusted).
// Returned sequences are in transcription (strand) order with array[0] being the first intron.
// The optional second offset argument is onlu used if the linear sequence is not a complete contig/chromosome.
std::vector<kgl::DNA5SequenceCoding>
kgl::DNA5SequenceLinear::intronArraySequence( const std::shared_ptr<const TranscriptionSequence>&,
                                              ContigOffset_t) const {
  return std::vector<kgl::DNA5SequenceCoding>();
//  return SequenceOffset::refIntronArraySequence(coding_seq_ptr, *this, 0, 0, contig_offset);

}
