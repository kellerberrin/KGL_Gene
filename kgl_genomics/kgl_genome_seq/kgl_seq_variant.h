//
// Created by kellerberrin on 3/01/18.
//

#ifndef KGL_SEQ_VARIANT_H
#define KGL_SEQ_VARIANT_H


#include "kgl_variant_db.h"
#include "kgl_mutation_db.h"

#include <map>
#include <memory>
#include <vector>
#include <sstream>



namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variant Mutation Functionality.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

private:


};


}   // end namespace


#endif //KGL_SEQ_VARIANT_H
