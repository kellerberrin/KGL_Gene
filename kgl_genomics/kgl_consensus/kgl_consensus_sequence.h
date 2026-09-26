//
// Module 3: consensus sequence construction.
//

#ifndef KGL_CONSENSUS_SEQUENCE_H
#define KGL_CONSENSUS_SEQUENCE_H

#include "kgl_consensus_offset.h"
#include "kgl_sequence.h"

#include <expected>
#include <optional>

namespace kellerberrin::genome {

struct ConsensusOptions {

  bool verify_reference{true};   // consumed bases must equal the variant reference

};

// The modified sequence and its mapping back to the reference.
//
//   bases()            : n bases, local coordinates [0, n).
//   origin()           : the supplied reference interval [R0, R1).
//   retainedReference(): [W0, R1); W0 > R0 only when an upstream delete clips the 5' end.
//   modifiedInterval() : [W0, W0 + n); always the same size as bases().
//                        Upstream inserts and SNPs have no effect on this interval.
class ConsensusSequence {

public:

  ConsensusSequence(ConsensusSequence&&) noexcept = default;
  ConsensusSequence& operator=(ConsensusSequence&&) noexcept = default;
  ConsensusSequence(const ConsensusSequence&) = delete;
  ConsensusSequence& operator=(const ConsensusSequence&) = delete;
  ~ConsensusSequence() = default;

  [[nodiscard]] const DNA5SequenceLinear& bases() const noexcept { return bases_; }
  [[nodiscard]] const OpenRightUnsigned& origin() const noexcept { return origin_; }
  [[nodiscard]] const OpenRightUnsigned& retainedReference() const noexcept { return retained_reference_; }
  [[nodiscard]] const OpenRightUnsigned& modifiedInterval() const noexcept { return modified_interval_; }
  [[nodiscard]] SignedOffset_t lengthDelta() const noexcept { return length_delta_; }
  [[nodiscard]] bool empty() const noexcept { return bases_.empty(); }

  // Moves the modified bases out (compatibility layer); the object is left empty.
  [[nodiscard]] DNA5SequenceLinear takeBases() noexcept { return std::move(bases_); }

  // Reference interval -> consensus sub-sequence (exon and buffer slicing).
  // nullopt when the interval lies outside the retained reference window;
  // an empty sequence when the interval was wholly deleted.
  [[nodiscard]] std::optional<DNA5SequenceLinear> slice(const OpenRightUnsigned& reference_interval) const;
  [[nodiscard]] bool isDeleted(const OpenRightUnsigned& reference_interval) const;

private:

  friend std::expected<ConsensusSequence, ConsensusError>
  buildConsensus(const DNA5SequenceLinear& reference,
                 const OpenRightUnsigned& origin,
                 const OffsetAccounting& accounting,
                 const ConsensusOptions& options);

  ConsensusSequence() = default;

  DNA5SequenceLinear bases_;
  OpenRightUnsigned origin_{0, 0};
  OpenRightUnsigned retained_reference_{0, 0};
  OpenRightUnsigned modified_interval_{0, 0};
  SignedOffset_t length_delta_{0};
  std::vector<Edit> applied_;   // edits that intersect the origin window

};

// `reference` holds the bases of `origin`: base i is contig offset origin.lower() + i.
// Applies the accounting's edits that can affect the window. Upstream inserts and SNPs
// are ignored; upstream deletes clip the retained reference start W0.
[[nodiscard]] std::expected<ConsensusSequence, ConsensusError>
buildConsensus(const DNA5SequenceLinear& reference,
               const OpenRightUnsigned& origin,
               const OffsetAccounting& accounting,
               const ConsensusOptions& options = {});

// Whole sequence convention (chromosome or scaffold): origin = [0, sequence.length()).
[[nodiscard]] std::expected<ConsensusSequence, ConsensusError>
buildConsensus(const DNA5SequenceLinear& sequence,
               const OffsetAccounting& accounting,
               const ConsensusOptions& options = {});

} // namespace

#endif //KGL_CONSENSUS_SEQUENCE_H
