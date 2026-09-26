//
// Module 2: interval mapping and offset accounting.
// Pure interval mathematics: no sequence type is used or included.
//

#ifndef KGL_CONSENSUS_OFFSET_H
#define KGL_CONSENSUS_OFFSET_H

#include "kgl_consensus_select.h"
#include "kgl_variant_db.h"

#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <vector>

namespace kellerberrin::genome {

// The kind of a resolved edit.
enum class EditKind : std::uint8_t { Snp, Insert, Delete };

// How overlapping or adjacent deletes in a resolved edit list are handled.
enum class DeleteOverlapPolicy : std::uint8_t { Union, Drop, Error };

// One resolved edit in contig coordinates.
//   SNP    : consumes [offset, offset+1), payload length 1.
//   Delete : consumes [offset+1, offset+1+n), payload empty.
//   Insert : consumes nothing (begin == end == offset+1), payload length n.
struct Edit {

  EditKind kind{EditKind::Snp};
  ContigOffset_t begin{0};
  ContigOffset_t end{0};
  std::shared_ptr<const Variant> variant;

  [[nodiscard]] ContigSize_t consumed() const noexcept { return end - begin; }
  [[nodiscard]] ContigSize_t inserted() const noexcept;
  // Alternate bases to emit: the alternate base for a SNP, the inserted bases
  // (alternate minus the shared anchor) for an insert, empty for a delete.
  // The variant is immutable so the view is stable for the edit's lifetime.
  [[nodiscard]] std::string_view insertedPayload() const noexcept;

};

// A globally resolved, canonical edit schedule for one contig.
// Post-conditions: edits are sorted by begin, disjoint, and shadow-pruned.
struct AppliedVariants {

  ContigId_t contig_id;
  std::vector<Edit> edits;
  FilteredVariantStats stats;

};

// Resolve a selected variant collection into a global edit schedule.
// Deletes are made disjoint first (per the overlap policy); non-delete edits that
// fall inside a delete are pruned (the "shadow of a delete" rule).
[[nodiscard]] std::expected<AppliedVariants, ConsensusError>
resolveVariants(const SelectedVariants& selected,
                DeleteOverlapPolicy overlap = DeleteOverlapPolicy::Union);

// Total, monotone mapping from reference contig coordinates to modified coordinates.
// An offset inside a deleted span collapses to the deletion boundary. The map is
// immutable once built and can be shared across threads.
class OffsetAccounting {

public:

  [[nodiscard]] static std::expected<OffsetAccounting, ConsensusError>
  build(const AppliedVariants& applied);

  // Reference -> modified. Total and monotone non-decreasing.
  [[nodiscard]] ContigOffset_t toModified(ContigOffset_t reference_offset) const noexcept;
  [[nodiscard]] OpenRightUnsigned toModified(const OpenRightUnsigned& reference_interval) const noexcept;

  // Modified position minus reference position. May be negative.
  [[nodiscard]] SignedOffset_t shift(ContigOffset_t reference_offset) const noexcept;

  [[nodiscard]] bool isDeleted(ContigOffset_t reference_offset) const noexcept;

  // Modified -> reference; nullopt inside an inserted payload.
  // At a deletion boundary the first retained reference base is returned.
  [[nodiscard]] std::optional<ContigOffset_t> toReference(ContigOffset_t modified_offset) const noexcept;

  [[nodiscard]] SignedOffset_t totalAdjust() const noexcept { return total_adjust_; }
  [[nodiscard]] std::span<const Edit> edits() const noexcept { return edits_; }
  [[nodiscard]] bool empty() const noexcept { return edits_.empty(); }

private:

  std::vector<Edit> edits_;                       // sorted by begin, disjoint
  std::vector<SignedOffset_t> shift_before_;      // accumulated delta before each edit
  std::vector<ContigOffset_t> modified_begin_;    // toModified(edits_[i].begin)
  SignedOffset_t total_adjust_{0};

};

// Map an arbitrary interval vector in one call.
// The image of [a, b) is [toModified(a), toModified(b)): unchanged before all edits,
// shifted after upstream edits, expanded over inserts, contracted or emptied by deletes.
[[nodiscard]] std::vector<OpenRightUnsigned>
modifyIntervals(const OffsetAccounting& accounting,
                std::span<const OpenRightUnsigned> reference_intervals);

} // namespace

#endif //KGL_CONSENSUS_OFFSET_H
