//
// Consensus library (glm-refactor polish of deepseek-refactor).
//
// Compatibility layer implementation.
//

#include "kgl_consensus_compat.h"
#include "kel_exec_env.h"

#include <algorithm>

namespace kgl = kellerberrin::genome;

bool kgl::AdjustedSequence::updateSequence(const std::shared_ptr<const ContigReference>& contig_ref_ptr,
                                           const SequenceVariantFilter& filtered_variants) {

  clear();

  if (not contig_ref_ptr) {

    return false;

  }

  const OpenRightUnsigned window = filtered_variants.sequenceInterval();
  region_ = window;

  // Build a SelectedVariants directly from the facade's filtered map: the facade
  // performs the window and policy filtering; the core resolves and applies.
  SelectedVariants selected;
  selected.stats = filtered_variants.filterStatistics();

  auto variant_contig = std::make_shared<ContigDB>(contig_ref_ptr->contigId());
  for (auto const& [offset, variant_ptr] : filtered_variants.offsetVariantMap()) {

    if (not variant_contig->addVariant(variant_ptr)) { return false; }

  }
  selected.variants = std::move(variant_contig);

  auto applied_opt = resolveVariants(selected);
  if (not applied_opt) {

    ExecEnv::log().warn("AdjustedSequence; resolution failed: {}", toString(applied_opt.error()));
    return false;

  }

  auto accounting_opt = OffsetAccounting::build(*applied_opt);
  if (not accounting_opt) {

    ExecEnv::log().warn("AdjustedSequence; accounting failed: {}", toString(accounting_opt.error()));
    return false;

  }

  auto reference_opt = contig_ref_ptr->sequence().subSequence(window);
  if (not reference_opt) {

    ExecEnv::log().warn("AdjustedSequence; cannot extract reference window: {}", window.toString());
    return false;

  }
  original_ = std::move(reference_opt.value());

  auto consensus_opt = buildConsensus(original_, window, *accounting_opt);
  if (not consensus_opt) {

    ExecEnv::log().warn("AdjustedSequence; consensus failed: {}", toString(consensus_opt.error()));
    return false;

  }
  consensus_.emplace(std::move(consensus_opt.value()));

  return true;

}

kgl::DualSeqOpt kgl::AdjustedSequence::moveSequenceClear() {

  if (not consensus_) {

    ExecEnv::log().warn("AdjustedSequence; no valid modified sequence available");
    return std::nullopt;

  }

  std::pair<DNA5SequenceLinear, DNA5SequenceLinear> pair{std::move(original_), consensus_->takeBases()};
  clear();
  return pair;

}

void kgl::AdjustedSequence::clear() {

  consensus_.reset();
  original_.clear();
  region_ = {0, 0};

}

std::optional<kgl::DNA5SequenceLinear>
kgl::AdjustedSequence::modifiedSubSequence(const OpenRightUnsigned& sub_interval) const {

  if (not consensus_) {

    ExecEnv::log().warn("AdjustedSequence; no valid modified sequence available");
    return std::nullopt;

  }

  return consensus_->slice(sub_interval);

}

std::optional<kgl::DNA5SequenceLinear>
kgl::AdjustedSequence::originalSubSequence(const OpenRightUnsigned& sub_interval) const {

  if (not consensus_) {

    ExecEnv::log().warn("AdjustedSequence; no valid reference sequence available");
    return std::nullopt;

  }

  if (not region_.containsInterval(sub_interval)) {

    ExecEnv::log().warn("AdjustedSequence; sub interval: {} is not contained in region: {}",
                        sub_interval.toString(), region_.toString());
    return std::nullopt;

  }

  return original_.subSequence({sub_interval.lower() - region_.lower(),
                                sub_interval.upper() - region_.lower()});

}
