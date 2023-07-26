//
// Created by kellerberrin on 26/07/23.
//

#ifndef KGL_MUTATION_OFFSET_H
#define KGL_MUTATION_OFFSET_H

#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"
#include "kgl_mutation_variant.h"
#include "kgl_mutation_db.h"


namespace kellerberrin::genome {   //  organization::project level namespace



class MutationOffset {

public:

  MutationOffset() = delete;
  ~MutationOffset() = delete;

  [[nodiscard]] static bool getSortedVariants( const std::shared_ptr<const GenomeDB>& genome_ptr,
                                               ContigId_t contig_id,
                                               VariantPhase phase,
                                               ContigOffset_t start,
                                               ContigOffset_t end,
                                               OffsetVariantMap &variant_map);

  [[nodiscard]] static bool getSortedVariants( const std::shared_ptr<const ContigDB>& contig_ptr,
                                               VariantPhase phase,
                                               ContigOffset_t start,
                                               ContigOffset_t end,
                                               OffsetVariantMap& variant_map);

private:

  void static checkUpstreamDeletion(OffsetVariantMap& variant_map);

};





} // Namespace

#endif //KGL_MUTATION_OFFSET_H
