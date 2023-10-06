//
// Created by kellerberrin on 26/07/23.
//

#ifndef KGL_SEQ_OFFSET_H
#define KGL_SEQ_OFFSET_H

#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"
#include "kgl_seq_variant.h"
#include "kgl_mutation_db.h"


namespace kellerberrin::genome {   //  organization::project level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Static functions return canonical variants over a specified region [start, end) ready to modify a DNA sequence
// over the same region and contig_ref_ptr.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SequenceVariantFilter {

public:

  SequenceVariantFilter(const std::shared_ptr<const ContigDB>& contig_ptr,
                        const OpenRightUnsigned& sequence_interval) {

    canonicalVariants( contig_ptr, sequence_interval);

  };
  ~SequenceVariantFilter() = default;

  // Returns canonical variants from the specified contig_ref_ptr and region specified as [start, end).
  // Note that delete indels can be have offset() < start.
  [[nodiscard]] static bool getSortedVariants( const std::shared_ptr<const GenomeDB>& genome_ptr,
                                               ContigId_t contig_id,
                                               VariantPhase phase,
                                               ContigOffset_t start,
                                               ContigOffset_t end,
                                               OffsetVariantMap &variant_map);



  [[nodiscard]] const OffsetVariantMap& offsetVariantMap() const { return offset_variant_map_; }
  [[nodiscard]] size_t duplicateVariants() const { return duplicate_variants_; }
  [[nodiscard]] size_t downstreamDelete() const { return downstream_delete_; }
  [[nodiscard]] size_t sequenceVariants() const { return offset_variant_map_.size() + duplicate_variants_ + downstream_delete_; }

private:

  OffsetVariantMap offset_variant_map_;
  size_t duplicate_variants_{0};
  size_t downstream_delete_{0};
  // A margin to account for the change in offsets when converting to canonical variants.
  // Canonical offsets always increase.
  constexpr static const SignedOffset_t NUCLEOTIDE_CANONICAL_MARGIN{200} ;

  // Returns a map of unique canonical variants
  // Also returns the number of multiple variants found at each offset which are filtered to a single variant.
  void canonicalVariants( const std::shared_ptr<const ContigDB>& contig_ptr, const OpenRightUnsigned& sequence_interval);
  // Returns a map of unique canonical variants
  // Also returns the number of multiple variants found at each offset which are filtered to a single variant.
  [[nodiscard]] static std::tuple<OffsetVariantMap, size_t, size_t> getCanonicalVariants( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                                          const OpenRightUnsigned& sequence_interval);


};





} // Namespace

#endif //KGL_SEQ_OFFSET_H
