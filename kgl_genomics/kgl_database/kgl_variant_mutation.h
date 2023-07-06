//
// Created by kellerberrin on 3/01/18.
//

#ifndef KGL_VARIANT_MUTATION_H
#define KGL_VARIANT_MUTATION_H

//
// Created by kellerberrin on 31/10/17.
//

#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_variant.h"
#include "kgl_variant_mutation_offset.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variant Mutation Functionality - Implements mutation functionality
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetVariantMap = std::map<ContigOffset_t, std::shared_ptr<const Variant>>;


class VariantMutation {

public:

  explicit VariantMutation() = default;
  ~VariantMutation() = default;

  // The coding variants in the variant_map are used to mutate the dna_sequence.
  [[nodiscard]] bool mutateDNA( const OffsetVariantMap &variant_map,
                                const std::shared_ptr<const ContigReference>& contig_ptr,
                                const std::shared_ptr<const TranscriptionSequence>& coding_sequence_ptr,
                                DNA5SequenceCoding &dna_sequence);

  [[nodiscard]] bool mutateDNA( const OffsetVariantMap &insert_variant_map,
                                const std::shared_ptr<const ContigReference>& contig_ptr,
                                ContigOffset_t contig_offset,
                                ContigSize_t sequence_size,
                                DNA5SequenceLinear& dna_sequence_ptr);

  [[nodiscard]] bool mutateDNA( const OffsetVariantMap &variant_map,
                                const std::shared_ptr<const ContigReference>& contig_ptr,
                                DNA5SequenceContig& contig_sequence_ptr);

private:

  VariantMutationOffset variant_mutation_offset_;

  [[nodiscard]] static bool mutateSequence( const DNA5SequenceLinear& reference,
                                            const DNA5SequenceLinear& alternate,
                                            ContigOffset_t canonical_offset,
                                            SignedOffset_t offset_adjust,
                                            DNA5SequenceLinear& dna_sequence,
                                            SignedOffset_t& sequence_size_modify);

  // Mutate a sequence by adding and subtracting subsequences at the designated canonical_offset
  [[nodiscard]] static bool performMutation( ContigOffset_t canonical_offset,
                                             DNA5SequenceLinear& mutated_sequence,
                                             const DNA5SequenceLinear& delete_subsequence,
                                             const DNA5SequenceLinear& add_subsequence,
                                             SignedOffset_t& sequence_size_modify);

  [[nodiscard]] static bool preceedingMutation( const DNA5SequenceLinear& reference,
                                                const DNA5SequenceLinear& alternate,
                                                SignedOffset_t adjusted_offset,
                                                DNA5SequenceLinear& dna_sequence,
                                                SignedOffset_t& sequence_size_modify);

};


}   // end namespace


#endif //READSAMFILE_KGL_VARIANT_MUTATION_H
