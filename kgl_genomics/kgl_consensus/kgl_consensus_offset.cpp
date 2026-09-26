//
// Module 2 implementation: interval mapping and offset accounting.
//

#include "kgl_consensus_offset.h"

#include "kel_exec_env.h"

#include <algorithm>
#include <ranges>

namespace kgl = kellerberrin::genome;

std::string_view kgl::toString(ConsensusError error) noexcept {

  switch (error) {

    case ConsensusError::NullContig: return "null contig";
    case ConsensusError::EmptyReference: return "empty reference";
    case ConsensusError::ReferenceOriginMismatch: return "reference/origin size mismatch";
    case ConsensusError::VariantNotCanonical: return "variant not canonical";
    case ConsensusError::VariantDuplicateLocus: return "duplicate variant locus";
    case ConsensusError::DeletesOverlap: return "overlapping deletes";
    case ConsensusError::ReferenceMismatch: return "reference base mismatch";
    case ConsensusError::LengthInvariant: return "length invariant violated";

  }

  return "unknown consensus error";

}

////////////////////////////////////////////////////////////////////////////////
// Edit helpers
////////////////////////////////////////////////////////////////////////////////

kgl::ContigSize_t kgl::Edit::inserted() const noexcept {

  if (kind == EditKind::Insert) {

    return variant->alternateSize() - variant->referenceSize();

  }
  if (kind == EditKind::Snp) {

    return 1;

  }
  return 0;

}

std::string_view kgl::Edit::insertedPayload() const noexcept {

  if (kind == EditKind::Snp) {

    return variant->alternate().getStringView().substr(0, 1);

  }
  if (kind == EditKind::Insert) {

    const ContigSize_t anchor = variant->referenceSize();
    return variant->alternate().getStringView().substr(anchor);

  }
  return {};

}

////////////////////////////////////////////////////////////////////////////////
// resolveVariants
////////////////////////////////////////////////////////////////////////////////

namespace {

using namespace kellerberrin;

// The reference span a variant claims for shadow purposes.
//   SNP    : [o, o+1)
//   Delete : [d0, d1)
//   Insert : [i_0-1, i_0+1) (anchor and insertion point)
kellerberrin::OpenRightUnsigned memberSpan(const genome::Edit& edit) {

  if (edit.kind == genome::EditKind::Insert) {

    return {edit.begin - 1, edit.begin + 1};

  }
  return {edit.begin, edit.end};

}

genome::Edit makeEdit(const std::shared_ptr<const genome::Variant>& variant_ptr) {

  const auto [variant_type, modify_interval] = variant_ptr->modifyInterval();

  genome::Edit edit;
  edit.variant = variant_ptr;

  switch (variant_type) {

    case genome::VariantType::SNP:
      edit.kind = genome::EditKind::Snp;
      edit.begin = variant_ptr->offset();
      edit.end = edit.begin + 1;
      break;

    case genome::VariantType::INDEL_DELETE:
      edit.kind = genome::EditKind::Delete;
      edit.begin = modify_interval.lower();
      edit.end = modify_interval.upper();
      break;

    case genome::VariantType::INDEL_INSERT:
      edit.kind = genome::EditKind::Insert;
      edit.begin = modify_interval.lower();
      edit.end = edit.begin;
      break;

  }

  return edit;

}

} // namespace

std::expected<kgl::AppliedVariants, kgl::ConsensusError>
kgl::resolveVariants(const SelectedVariants& selected, DeleteOverlapPolicy overlap) {

  if (not selected.variants) {

    return std::unexpected(ConsensusError::NullContig);

  }

  AppliedVariants applied;
  applied.contig_id = selected.variants->contigId();
  applied.stats = selected.stats;

  std::vector<Edit> edits;
  bool failed{false};
  selected.variants->processAll([&edits, &failed](const std::shared_ptr<const Variant>& variant_ptr) {

    if (not variant_ptr->isCanonical()) {

      kellerberrin::ExecEnv::log().error("consensus::resolveVariants; variant NOT canonical: {}", variant_ptr->HGVS());
      failed = true;
      return false;

    }
    edits.push_back(makeEdit(variant_ptr));
    return true;

  });

  if (failed) {

    return std::unexpected(ConsensusError::VariantNotCanonical);

  }

  std::ranges::sort(edits, [](const Edit& lhs, const Edit& rhs) { return lhs.begin < rhs.begin; });

  // Partition deletes; resolve them into a disjoint set first.
  std::vector<Edit> deletes;
  std::vector<Edit> others;
  for (auto const& edit : edits) {

    (edit.kind == EditKind::Delete ? deletes : others).push_back(edit);

  }

  std::vector<Edit> disjoint_deletes;
  for (auto const& del : deletes) {

    if (disjoint_deletes.empty() or del.begin > disjoint_deletes.back().end) {

      disjoint_deletes.push_back(del);
      continue;

    }

    switch (overlap) {

      case DeleteOverlapPolicy::Union: {
        Edit& previous = disjoint_deletes.back();
        if (del.end > previous.end) {

          previous.end = del.end;

        }
        previous.variant = del.variant;
      }
      break;

      case DeleteOverlapPolicy::Drop:
        break;

      case DeleteOverlapPolicy::Error:
        kellerberrin::ExecEnv::log().warn("consensus::resolveVariants; overlapping deletes at {} and {}",
                            disjoint_deletes.back().begin, del.begin);
        return std::unexpected(ConsensusError::DeletesOverlap);

    }

  }

  // Shadow pruning: a non-delete edit whose member span intersects a delete is erased.
  std::vector<Edit> surviving;
  surviving.reserve(others.size());
  for (auto const& edit : others) {

    const kellerberrin::OpenRightUnsigned member = memberSpan(edit);
    bool shadowed{false};
    for (auto const& del : disjoint_deletes) {

      if (OpenRightUnsigned{del.begin, del.end}.intersects(member)) {

        shadowed = true;
        break;

      }

    }

    if (shadowed) {

      ExecEnv::log().info("consensus::resolveVariants; edit at offset: {} shadowed by a delete, variant: {}",
                          edit.begin, edit.variant->HGVS());
      ++applied.stats.upstream_deleted_;

    } else {

      surviving.push_back(edit);

    }

  }

  surviving.insert(surviving.end(), disjoint_deletes.begin(), disjoint_deletes.end());
  std::ranges::sort(surviving, [](const Edit& lhs, const Edit& rhs) { return lhs.begin < rhs.begin; });

  // One edit per locus after pruning. (Competing alleles were resolved during selection.)
  for (size_t index = 1; index < surviving.size(); ++index) {

    if (surviving[index].begin == surviving[index - 1].begin) {

      kellerberrin::ExecEnv::log().warn("consensus::resolveVariants; duplicate locus at offset: {}", surviving[index].begin);
      return std::unexpected(ConsensusError::VariantDuplicateLocus);

    }

  }

  applied.edits = std::move(surviving);

  return applied;

}

////////////////////////////////////////////////////////////////////////////////
// OffsetAccounting
////////////////////////////////////////////////////////////////////////////////

std::expected<kgl::OffsetAccounting, kgl::ConsensusError>
kgl::OffsetAccounting::build(const AppliedVariants& applied) {

  OffsetAccounting accounting;
  accounting.edits_ = applied.edits;

  for (size_t index = 1; index < accounting.edits_.size(); ++index) {

    if (accounting.edits_[index].begin < accounting.edits_[index - 1].end) {

      return std::unexpected(ConsensusError::DeletesOverlap);

    }

  }

  SignedOffset_t cumulative{0};
  accounting.shift_before_.reserve(accounting.edits_.size());
  accounting.modified_begin_.reserve(accounting.edits_.size());
  for (auto const& edit : accounting.edits_) {

    accounting.shift_before_.push_back(cumulative);
    accounting.modified_begin_.push_back(
        static_cast<ContigOffset_t>(static_cast<SignedOffset_t>(edit.begin) + cumulative));
    cumulative += static_cast<SignedOffset_t>(edit.inserted()) - static_cast<SignedOffset_t>(edit.consumed());

  }

  accounting.total_adjust_ = cumulative;

  return accounting;

}

namespace {

// Index of the first edit beginning after the offset. Both edit lists are sorted
// by begin, so this is also the number of edits with begin <= offset (O(log V)).
[[nodiscard]] size_t indexOfFirstEditAfter(const std::vector<genome::Edit>& edits, genome::ContigOffset_t offset) {

  const auto first_after = std::ranges::upper_bound(edits, offset, std::ranges::less{}, &genome::Edit::begin);
  return static_cast<size_t>(std::ranges::distance(edits.begin(), first_after));

}

} // namespace

kgl::ContigOffset_t kgl::OffsetAccounting::toModified(ContigOffset_t reference_offset) const noexcept {

  const size_t count = indexOfFirstEditAfter(edits_, reference_offset);
  if (count == 0) {

    return reference_offset;

  }

  const size_t index = count - 1;
  const Edit& edit = edits_[index];

  if (edit.kind == EditKind::Delete and reference_offset < edit.end) {

    return modified_begin_[index];   // collapse to the deletion boundary

  }

  // Insert boundary convention (matches ConsensusSequence::slice and the reference
  // implementation): an offset exactly at the insertion point maps to the payload
  // START, so the payload is included only by intervals whose lower endpoint is the
  // insertion point, never by intervals ending there. Without this, an interval
  // ending at the point would absorb the payload and disagree with slice().
  if (edit.kind == EditKind::Insert and reference_offset == edit.begin) {

    return static_cast<ContigOffset_t>(static_cast<SignedOffset_t>(reference_offset) + shift_before_[index]);

  }

  const SignedOffset_t cumulative =
      shift_before_[index] + static_cast<SignedOffset_t>(edit.inserted()) - static_cast<SignedOffset_t>(edit.consumed());
  return static_cast<ContigOffset_t>(static_cast<SignedOffset_t>(reference_offset) + cumulative);

}

kellerberrin::OpenRightUnsigned kgl::OffsetAccounting::toModified(const kellerberrin::OpenRightUnsigned& reference_interval) const noexcept {

  return {toModified(reference_interval.lower()), toModified(reference_interval.upper())};

}

kgl::SignedOffset_t kgl::OffsetAccounting::shift(ContigOffset_t reference_offset) const noexcept {

  return static_cast<SignedOffset_t>(toModified(reference_offset)) - static_cast<SignedOffset_t>(reference_offset);

}

bool kgl::OffsetAccounting::isDeleted(ContigOffset_t reference_offset) const noexcept {

  const size_t count = indexOfFirstEditAfter(edits_, reference_offset);
  if (count == 0) {

    return false;

  }

  const Edit& edit = edits_[count - 1];
  return edit.kind == EditKind::Delete and reference_offset < edit.end;

}

std::optional<kgl::ContigOffset_t> kgl::OffsetAccounting::toReference(ContigOffset_t modified_offset) const noexcept {

  // Index of the last edit whose modified begin is <= modified_offset.
  const auto first_after = std::ranges::upper_bound(modified_begin_, modified_offset);
  if (first_after == modified_begin_.begin()) {

    return modified_offset;

  }

  const size_t index = static_cast<size_t>(std::ranges::distance(modified_begin_.begin(), first_after)) - 1;
  const Edit& edit = edits_[index];

  if (edit.kind == EditKind::Insert and modified_offset < modified_begin_[index] + edit.inserted()) {

    return std::nullopt;   // no reference base inside an inserted payload

  }

  const SignedOffset_t cumulative =
      shift_before_[index] + static_cast<SignedOffset_t>(edit.inserted()) - static_cast<SignedOffset_t>(edit.consumed());
  return static_cast<ContigOffset_t>(static_cast<SignedOffset_t>(modified_offset) - cumulative);

}

std::vector<kellerberrin::OpenRightUnsigned>
kgl::modifyIntervals(const OffsetAccounting& accounting,
                     std::span<const OpenRightUnsigned> reference_intervals) {

  std::vector<OpenRightUnsigned> modified;
  modified.reserve(reference_intervals.size());
  for (auto const& interval : reference_intervals) {

    modified.push_back(accounting.toModified(interval));

  }

  return modified;

}
