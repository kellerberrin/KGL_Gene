//
// Created by kellerberrin on 23/09/23.
//

#include "kgl_genome_seq/kgl_seq_variant.h"


namespace kgl = kellerberrin::genome;

// Delete ASAP

// The coding variants in the variant_map are used to mutate the dna_sequence.
 bool kgl::VariantMutation::mutateDNA( const OffsetVariantMap&,
                              const std::shared_ptr<const ContigReference>&,
                              const std::shared_ptr<const TranscriptionSequence>&,
                              DNA5SequenceCoding&) { return true; }

bool kgl::VariantMutation::mutateDNA( const OffsetVariantMap &,
                              const std::shared_ptr<const ContigReference>&,
                              ContigOffset_t,
                              ContigSize_t,
                              DNA5SequenceLinear&) { return true; }
