//
// Created by kellerberrin on 3/01/18.
//

#ifndef KGL_VARIANT_MUTATION_H
#define KGL_VARIANT_MUTATION_H

//
// Created by kellerberrin on 31/10/17.
//

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
                                DNA5SequenceCoding &dna_sequence) { return true; }

  [[nodiscard]] bool mutateDNA( const OffsetVariantMap &insert_variant_map,
                                const std::shared_ptr<const ContigReference>& contig_ptr,
                                ContigOffset_t contig_offset,
                                ContigSize_t sequence_size,
                                DNA5SequenceLinear& dna_sequence_ptr) { return true; }

private:


};


}   // end namespace


#endif //READSAMFILE_KGL_VARIANT_MUTATION_H
