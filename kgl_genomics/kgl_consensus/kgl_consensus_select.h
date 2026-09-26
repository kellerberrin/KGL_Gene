//
// Module 1: variant selection and filtering.
// Produces a canonical, policy-filtered, one-variant-per-locus ContigDB.
//

#ifndef KGL_CONSENSUS_SELECT_H
#define KGL_CONSENSUS_SELECT_H

#include "kgl_variant_db_contig.h"
#include "kel_interval_unsigned.h"

#include <cstdint>
#include <expected>
#include <map>
#include <memory>
#include <string_view>

namespace kellerberrin::genome {

// Errors returned by the consensus modules.
enum class ConsensusError : std::uint8_t {
  NullContig,
  EmptyReference,
  ReferenceOriginMismatch,
  VariantNotCanonical,
  VariantDuplicateLocus,
  DeletesOverlap,
  ReferenceMismatch,
  LengthInvariant,
};

[[nodiscard]] std::string_view toString(ConsensusError error) noexcept;

// Statistics for the variants selected to modify a region.
// Field names and types are unchanged from the reference implementation; the struct
// is embedded by value in a public client statistics object.
struct FilteredVariantStats {

  size_t total_interval_variants_{0};
  size_t total_snp_variants_{0};
  size_t total_frame_shift_{0};
  size_t non_unique_count_{0};
  size_t upstream_deleted_{0};

};

using OffsetVariantMap = std::map<ContigOffset_t, std::shared_ptr<const Variant>>;

// The type of sequence variant filter selected by the client.
enum class SeqVariantFilterType { DEFAULT_SEQ_FILTER, HIGHEST_FREQ_VARIANT, FRAMESHIFT_ADJUSTED, SNP_ADJUSTED };

// A canonical, policy-filtered, one-variant-per-locus variant collection for a whole contig.
// Shadow pruning and window geometry are the responsibility of the offset module, which
// sees the whole contig and therefore requires no upstream margin heuristic.
struct SelectedVariants {

  std::shared_ptr<const ContigDB> variants;
  FilteredVariantStats stats;

};

// Canonicalise, apply the policy and resolve multiple variants per locus.
// Contig-global selection; not used by the transcript client (which is window-scoped)
// but retained for whole-contig / whole-chromosome consensus builds.
// Non-canonical variants are logged and skipped; they are never silently applied.
[[nodiscard]] std::expected<SelectedVariants, ConsensusError>
selectVariants(const std::shared_ptr<const ContigDB>& raw_variants,
               SeqVariantFilterType policy = SeqVariantFilterType::DEFAULT_SEQ_FILTER);

// Windowed selection: the pipeline used by transcript clients.
// Reproduces the reference/constrained collection exactly:
//   margin -> ContigModifyFilter -> statistics -> policy -> per-locus resolution
//   -> shadow pruning by upstream deletes.
// The returned map holds the variants to apply; the statistics are window-scoped.
[[nodiscard]] std::expected<SelectedVariants, ConsensusError>
selectWindowVariants(const std::shared_ptr<const ContigDB>& raw_variants,
                     const OpenRightUnsigned& window,
                     SeqVariantFilterType policy = SeqVariantFilterType::DEFAULT_SEQ_FILTER);

// Compatibility facade preserving the historical public surface used by client code.
// Internally applies the same windowed pipeline as the reference implementation.
class SequenceVariantFilter {

public:

  SequenceVariantFilter(const std::shared_ptr<const ContigDB>& contig_ptr,
                        const OpenRightUnsigned& sequence_interval,
                        SeqVariantFilterType seq_filter_type = SeqVariantFilterType::DEFAULT_SEQ_FILTER);
  ~SequenceVariantFilter() = default;

  [[nodiscard]] const OpenRightUnsigned& sequenceInterval() const { return sequence_interval_; }
  [[nodiscard]] SeqVariantFilterType sequenceFilterType() const { return sequence_filter_type_; }
  [[nodiscard]] const OffsetVariantMap& offsetVariantMap() const { return offset_variant_map_; }
  [[nodiscard]] const FilteredVariantStats& filterStatistics() const { return variant_filter_stats_; }

private:

  const OpenRightUnsigned sequence_interval_;
  SeqVariantFilterType sequence_filter_type_;
  OffsetVariantMap offset_variant_map_;
  FilteredVariantStats variant_filter_stats_;

  // A margin to account for the change in offsets when converting to canonical variants.
  constexpr static const SignedOffset_t NUCLEOTIDE_CANONICAL_MARGIN_{200};

  void collect(const std::shared_ptr<const ContigDB>& contig_ptr, const OpenRightUnsigned& sequence_interval);

};

} // namespace

#endif //KGL_CONSENSUS_SELECT_H
