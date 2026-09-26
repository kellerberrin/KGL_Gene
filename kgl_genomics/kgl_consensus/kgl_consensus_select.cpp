//
// Module 1 implementation: variant selection and filtering.
//

#include "kgl_consensus_select.h"

#include "kgl_variant_filter_db_contig.h"
#include "kgl_variant_filter_db_offset.h"
#include "kgl_variant_filter_db_variant.h"
#include "kgl_variant_filter_coding.h"
#include "kel_exec_env.h"

#include <algorithm>
#include <set>

namespace kgl = kellerberrin::genome;

namespace {

// The apply key of a canonical variant in the historical offset map:
// an SNP is keyed at its offset, an indel at the offset the edit occurs (offset + 1).
[[nodiscard]] kgl::ContigOffset_t applyKey(const kgl::Variant& variant) {

  return variant.isSNP() ? variant.offset() : variant.offset() + 1;

}

// Canonicalise and reduce a raw contig to a canonical, one-variant-per-locus,
// policy-filtered ContigDB. The window is applied by the caller for the
// compatibility facade; selection itself is contig-global.
std::shared_ptr<const kgl::ContigDB>
canonicalContig(const std::shared_ptr<const kgl::ContigDB>& raw_ptr) {

  auto canonical_ptr = std::make_shared<kgl::ContigDB>(raw_ptr->contigId());
  raw_ptr->processAll([&canonical_ptr](const std::shared_ptr<const kgl::Variant>& variant_ptr) {

    if (not variant_ptr->isCanonical()) {

      // Convert to canonical form on demand and add.
      auto clone_ptr = variant_ptr->cloneCanonical();
      if (clone_ptr and clone_ptr->isCanonical()) {

        if (not canonical_ptr->addVariant(std::shared_ptr<const kgl::Variant>(std::move(clone_ptr)))) {

          kellerberrin::ExecEnv::log().error("consensus::selectVariants; could not add canonical variant: {}", variant_ptr->HGVS());

        }

      } else {

        kellerberrin::ExecEnv::log().warn("consensus::selectVariants; non-canonical variant skipped: {}", variant_ptr->HGVS());

      }

      return true;   // skip: a bad variant never aborts the whole contig scan

    }

    if (not canonical_ptr->addVariant(variant_ptr)) {

      kellerberrin::ExecEnv::log().error("consensus::selectVariants; could not add variant: {}", variant_ptr->HGVS());

    }

    return true;

  });

  return canonical_ptr;

}

} // namespace

std::expected<kgl::SelectedVariants, kgl::ConsensusError>
kgl::selectVariants(const std::shared_ptr<const ContigDB>& raw_variants, SeqVariantFilterType policy) {

  if (not raw_variants) {

    return std::unexpected(ConsensusError::NullContig);

  }

  SelectedVariants selected;
  selected.stats.total_interval_variants_ = raw_variants->variantCount();

  // Count SNPs and frameshifts before resolution (the frame-shift statistic is
  // policy-independent: it reports what was present, not what was applied).
  raw_variants->processAll([&selected](const std::shared_ptr<const Variant>& variant_ptr) {

    if (variant_ptr->isSNP()) {

      ++selected.stats.total_snp_variants_;

    } else if ((variant_ptr->modifyInterval().second.size() % 3) != 0) {

      ++selected.stats.total_frame_shift_;

    }
    return true;

  });

  auto canonical_ptr = canonicalContig(raw_variants);

  // Resolve multiple variants at each locus. HomozygousCodingFilter prefers
  // homozygous variants and falls back to frequency, as in the reference.
  auto unique_ptr = canonical_ptr->viewFilter(HomozygousCodingFilter());
  const size_t unique_count = unique_ptr->variantCount();

  // Apply the policy filter to the unique set.
  std::unique_ptr<ContigDB> filtered_ptr;
  switch (policy) {

    case SeqVariantFilterType::FRAMESHIFT_ADJUSTED:
      filtered_ptr = unique_ptr->viewFilter(NotFilter(FrameShiftFilter()));
      break;

    case SeqVariantFilterType::SNP_ADJUSTED:
      filtered_ptr = unique_ptr->viewFilter(SNPFilter());
      break;

    case SeqVariantFilterType::DEFAULT_SEQ_FILTER:
    case SeqVariantFilterType::HIGHEST_FREQ_VARIANT:
      filtered_ptr = std::move(unique_ptr);
      break;

  }

  // Non-unique count is the duplicate reduction only (policy removal is not duplication).
  selected.stats.non_unique_count_ = canonical_ptr->variantCount() - unique_count;
  selected.variants = std::shared_ptr<const ContigDB>(std::move(filtered_ptr));

  return selected;

}

////////////////////////////////////////////////////////////////////////////////
// Windowed selection (the transcript production path)
////////////////////////////////////////////////////////////////////////////////

namespace {

// Canonical-offset margin; mirrors the reference filter constant.
constexpr kellerberrin::genome::SignedOffset_t WINDOW_CANONICAL_MARGIN_{200};

} // namespace

std::expected<kgl::SelectedVariants, kgl::ConsensusError>
kgl::selectWindowVariants(const std::shared_ptr<const ContigDB>& raw_variants,
                          const OpenRightUnsigned& window,
                          SeqVariantFilterType policy) {

  if (not raw_variants) {

    return std::unexpected(ConsensusError::NullContig);

  }

  // Margin below the window for non-canonical conversion and upstream deletes,
  // exactly as the reference/constrained collection functions.
  const ContigOffset_t lower = std::max<SignedOffset_t>(
      0, static_cast<SignedOffset_t>(window.lower()) - WINDOW_CANONICAL_MARGIN_);

  auto region_ptr = raw_variants->viewFilter(ContigRegionFilter(lower, window.upper()));
  auto modify_ptr = region_ptr->viewFilter(ContigModifyFilter(window.lower(), window.upper()));

  // Counters before the policy filter (the intended statistics semantics: they
  // report what is present in the window, not what the policy retains).
  auto hetero_ptr = modify_ptr->viewFilter(HeterozygousFilter());
  FilteredVariantStats stats;
  stats.total_interval_variants_ = hetero_ptr->variantCount();
  stats.total_snp_variants_ = hetero_ptr->viewFilter(SNPFilter())->variantCount();
  stats.total_frame_shift_ = hetero_ptr->viewFilter(FrameShiftFilter())->variantCount();

  // Policy filter.
  std::unique_ptr<ContigDB> policy_ptr;
  switch (policy) {

    case SeqVariantFilterType::FRAMESHIFT_ADJUSTED:
      policy_ptr = modify_ptr->viewFilter(NotFilter(FrameShiftFilter()));
      break;

    case SeqVariantFilterType::SNP_ADJUSTED:
      policy_ptr = modify_ptr->viewFilter(SNPFilter());
      break;

    case SeqVariantFilterType::DEFAULT_SEQ_FILTER:
    case SeqVariantFilterType::HIGHEST_FREQ_VARIANT:
      policy_ptr = std::move(modify_ptr);
      break;

  }

  // Per-locus resolution, then shadow pruning by upstream deletes.
  auto unique_ptr = policy_ptr->viewFilter(HomozygousCodingFilter());
  auto surviving_ptr = unique_ptr->viewFilter(ContigUpstreamFilter());

  stats.upstream_deleted_ = unique_ptr->variantCount() - surviving_ptr->variantCount();

  // Non-unique count follows the reference: unique-unphased count of the policy set
  // minus the number of distinct apply keys in the surviving set.
  const size_t modify_count = policy_ptr->viewFilter(UniqueUnphasedFilter())->variantCount();
  std::set<ContigOffset_t> keys;
  for (auto const& [offset, offset_ptr] : surviving_ptr->getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      keys.insert(applyKey(*variant_ptr));

    }

  }
  stats.non_unique_count_ = modify_count > keys.size() ? modify_count - keys.size() : 0;

  SelectedVariants selected;
  selected.stats = stats;
  selected.variants = std::shared_ptr<const ContigDB>(std::move(surviving_ptr));

  return selected;

}

////////////////////////////////////////////////////////////////////////////////
// SequenceVariantFilter compatibility facade
////////////////////////////////////////////////////////////////////////////////

kgl::SequenceVariantFilter::SequenceVariantFilter(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                  const OpenRightUnsigned& sequence_interval,
                                                  SeqVariantFilterType seq_filter_type)
    : sequence_interval_(sequence_interval), sequence_filter_type_(seq_filter_type) {

  collect(contig_ptr, sequence_interval_);

}

void kgl::SequenceVariantFilter::collect(const std::shared_ptr<const ContigDB>& contig_ptr,
                                         const OpenRightUnsigned& sequence_interval) {

  if (not contig_ptr) {

    return;

  }

  auto selected_opt = selectWindowVariants(contig_ptr, sequence_interval, sequence_filter_type_);
  if (not selected_opt) {

    ExecEnv::log().warn("SequenceVariantFilter; selection failed: {}", toString(selected_opt.error()));
    return;

  }
  variant_filter_stats_ = selected_opt->stats;

  // Rebuild the historical offset-keyed map (SNP at offset, indel at offset + 1).
  for (auto const& [offset, offset_ptr] : selected_opt->variants->getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      offset_variant_map_.try_emplace(applyKey(*variant_ptr), variant_ptr);

    }

  }

}
