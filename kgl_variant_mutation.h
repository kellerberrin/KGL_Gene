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


  static void getMutationAlternatives(std::shared_ptr<const OffsetVariantMap> variant_map_ptr,
                                      std::vector<OffsetVariantMap> &variant_map_vector,
                                      size_t &alternative_count) {

    getMutationAlternatives(variant_map_ptr, variant_map_vector, alternative_count, MUTATION_SOFT_LIMIT_, MUTATION_HARD_LIMIT_);

  }


private:

  constexpr static size_t MUTATION_SOFT_LIMIT_ = 32;
  constexpr static size_t MUTATION_HARD_LIMIT_ = 128;

  VariantMutationOffset variant_mutation_offset_;

  // Returns a vector of alternative mutation paths based on the number of equal offset mutations in the coding variants.
  // There maybe more than one variant specified per offset.
  // If there are equal offset variants then we create alternative mutation paths.
  // This function is exponential. Alternative Mutations = 2 ^ (#equal offset variants).
  // A warning is issued if the soft limit is reached; default 32 alternatives (5 equal offset variants).
  // The number of variant paths is capped by the hard limit; default 128 alternatives (9 equal offset variants)
  static void getMutationAlternatives(std::shared_ptr<const OffsetVariantMap> variant_map_ptr,
                                      std::vector<OffsetVariantMap> &variant_map_vector,
                                      size_t &alternative_count,
                                      size_t soft_limit,
                                      size_t hard_limit);

  // Split the variant map into SNP, Delete and Insert Variants.
  static void SplitVariantMap(const OffsetVariantMap &variant_map,
                              OffsetVariantMap &snp_variant_map,
                              OffsetVariantMap &delete_variant_map,
                              OffsetVariantMap &insert_variant_map);

};


}   // namespace genome
}   // namespace kellerberrin


#endif //READSAMFILE_KGL_VARIANT_MUTATION_H
