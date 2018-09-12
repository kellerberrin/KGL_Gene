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
#include "kgl_variant_db.h"
#include "kgl_variant_mutation_offset.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variant Mutation Functionality - Implements mutation functionality
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantMutation {

public:

  explicit VariantMutation() = default;
  ~VariantMutation() = default;

  // The coding variants in the variant_map are used to mutate the dna_sequence.
  bool mutateDNA(const OffsetVariantMap &variant_map,
                 std::shared_ptr<const ContigFeatures> contig_ptr,
                 std::shared_ptr<const CodingSequence> coding_sequence_ptr,
                 std::shared_ptr<DNA5SequenceCoding> &dna_sequence_ptr);

  bool mutateDNA(const OffsetVariantMap &insert_variant_map,
                 ContigOffset_t contig_offset,
                 std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr);



private:

  VariantMutationOffset variant_mutation_offset_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //READSAMFILE_KGL_VARIANT_MUTATION_H
