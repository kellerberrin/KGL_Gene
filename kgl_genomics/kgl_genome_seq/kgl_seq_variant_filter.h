//
// Created by kellerberrin on 26/07/23.
//

#ifndef KGL_SEQ_VARIANT_FILTER_H
#define KGL_SEQ_VARIANT_FILTER_H

#include "kgl_variant_db_contig.h"
#include "kel_interval_unsigned.h"

namespace kellerberrin::genome {   //  organization::project level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Functions return canonical variants over a specified region [start, end) ready to modify a DNA sequence
// over the same region and contig_ref_ptr.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Object holds holds unique canonical variants for a specified region or transcript.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using OffsetVariantMap = std::map<ContigOffset_t, std::shared_ptr<const Variant>>;

// The type of sequence variant filter specified.
enum class SeqVariantFilterType { DEFAULT_SEQ_FILTER, HIGHEST_FREQ_VARIANT };

class SequenceVariantFilter {

public:

  SequenceVariantFilter(const std::shared_ptr<const ContigDB>& contig_ptr,
                        const OpenRightUnsigned& sequence_interval,
                        SeqVariantFilterType seq_filter_type = SeqVariantFilterType::DEFAULT_SEQ_FILTER)
    : sequence_interval_(sequence_interval), sequence_filter_type_(seq_filter_type) {

    selectFilterType( contig_ptr, sequence_interval_);

  };
  ~SequenceVariantFilter() = default;

  // The interval over which the variants are filtered. Note the there may be an upstream indel delete not in the interval.
  [[nodiscard]] const OpenRightUnsigned& sequenceInterval() const { return sequence_interval_; }
  // How the raw variants were filtered.
  [[nodiscard]] SeqVariantFilterType sequenceFilterType() const { return sequence_filter_type_; }
  // The filtered variants held in an offset keyed map.
  [[nodiscard]] const OffsetVariantMap& offsetVariantMap() const { return offset_variant_map_; }
  // Variants remove because they logically mutate the same section of DNA.
  [[nodiscard]] size_t duplicateVariants() const { return duplicate_variants_; }
  // Variants removed because they are downstream of an indel delete and modify deleted nucleotides.
  [[nodiscard]] size_t downstreamDelete() const { return downstream_delete_; }
  // All variants before filtering, includes any upstream indel delete that modifies the sequence interval.
  [[nodiscard]] size_t preFilterVariants() const { return offset_variant_map_.size() + duplicate_variants_ + downstream_delete_; }

private:

  const OpenRightUnsigned sequence_interval_;
  SeqVariantFilterType sequence_filter_type_;
  OffsetVariantMap offset_variant_map_;
  size_t duplicate_variants_{0};
  size_t downstream_delete_{0};

  // A margin to account for the change in offsets when converting to canonical variants.
  // Canonical offsets always increase.
  constexpr static const SignedOffset_t NUCLEOTIDE_CANONICAL_MARGIN{200} ;

  void selectFilterType(const std::shared_ptr<const ContigDB>& contig_ptr, const OpenRightUnsigned& sequence_interval);
  // Returns a map of unique canonical variants
  // Also returns the number of multiple variants found at each offset which are filtered to a single variant.
  void canonicalVariants( const std::shared_ptr<const ContigDB>& contig_ptr, const OpenRightUnsigned& sequence_interval);
  // Returns a map of unique canonical variants
  // Also returns the number of multiple variants found at each offset which are filtered to a single variant.
  [[nodiscard]] static std::tuple<OffsetVariantMap, size_t, size_t> getCanonicalVariants( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                                          const OpenRightUnsigned& sequence_interval);


};





} // Namespace

#endif //KGL_SEQ_VARIANT_FILTER_H
