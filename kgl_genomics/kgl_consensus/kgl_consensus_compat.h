//
// Consensus library (glm-refactor polish of deepseek-refactor).
//
// Compatibility layer: preserves the historical AdjustedSequence surface used by
// the two legacy clients. Implemented over modules 1-3.
//

#ifndef KGL_CONSENSUS_COMPAT_H
#define KGL_CONSENSUS_COMPAT_H

#include "kgl_consensus_select.h"
#include "kgl_consensus_sequence.h"
#include "kgl_genome_contig.h"

#include <optional>
#include <utility>

namespace kellerberrin::genome {

// Historical alias retained for client code.
using DualSeqOpt = std::optional<std::pair<DNA5SequenceLinear, DNA5SequenceLinear>>;

// Compatibility facade over the consensus core. The public surface matches the
// reference AdjustedSequence for the members the legacy clients use.
class AdjustedSequence {

public:

  AdjustedSequence() = default;
  ~AdjustedSequence() = default;

  // Update the region sequence using a pre-filtered variant set.
  [[nodiscard]] bool updateSequence(const std::shared_ptr<const ContigReference>& contig_ref_ptr,
                                    const SequenceVariantFilter& filtered_variants);

  // Move the reference and modified sequences and clear the object.
  [[nodiscard]] DualSeqOpt moveSequenceClear();

  [[nodiscard]] bool validModifiedSequence() const noexcept { return consensus_.has_value(); }

  // The zero-based modified sub-sequence for a contig interval (reference API).
  // Returns nullopt when the interval is outside the region; an empty sequence
  // when the interval was wholly deleted.
  [[nodiscard]] std::optional<DNA5SequenceLinear> modifiedSubSequence(const OpenRightUnsigned& sub_interval) const;

  // The zero-based unmodified sub-sequence for a contig interval (reference API).
  [[nodiscard]] std::optional<DNA5SequenceLinear> originalSubSequence(const OpenRightUnsigned& sub_interval) const;

  // The region interval supplied to updateSequence (reference API).
  [[nodiscard]] const OpenRightUnsigned& contigInterval() const noexcept { return region_; }

private:

  void clear();

  std::optional<ConsensusSequence> consensus_;
  DNA5SequenceLinear original_;
  OpenRightUnsigned region_{0, 0};

};

} // namespace

#endif //KGL_CONSENSUS_COMPAT_H
