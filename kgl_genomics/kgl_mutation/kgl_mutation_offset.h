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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Static functions return canonical variants over a specified region [start, end) ready to modify a DNA sequence
// over the same region and contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MutationOffset {

public:

  MutationOffset() = delete;
  ~MutationOffset() = delete;

  // Returns canonical variants from the specified contig and region specified as [start, end).
  // Note that delete indels can be have offset() < start.
  [[nodiscard]] static bool getSortedVariants( const std::shared_ptr<const GenomeDB>& genome_ptr,
                                               ContigId_t contig_id,
                                               VariantPhase phase,
                                               ContigOffset_t start,
                                               ContigOffset_t end,
                                               OffsetVariantMap &variant_map);

  [[nodiscard]] static std::pair<OffsetVariantMap, bool> getCanonicalVariants( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                               ContigOffset_t start,
                                                                               ContigOffset_t end);

private:

  // A margin to account for the change in offsets when converting to canonical variants.
  // Canonical offsets always increase.
  constexpr static const size_t NUCLEOTIDE_CANONICAL_MARGIN{200} ;

};





} // Namespace

#endif //KGL_MUTATION_OFFSET_H
