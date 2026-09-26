//
// Consensus library (glm-refactor polish of deepseek-refactor).
//
// Module 3 implementation: consensus sequence construction.
//

#include "kgl_consensus_sequence.h"
#include "kel_exec_env.h"

#include <algorithm>

namespace kgl = kellerberrin::genome;

namespace {

// True when a delete edit's span exactly matches its variant's canonical modify interval
// (i.e. it was not a union-merged delete, whose record no longer describes the full span).
bool deleteVariantMatches(const kgl::Edit& edit) {

  if (not edit.variant or edit.variant->variantType() != kgl::VariantType::INDEL_DELETE) {

    return false;

  }

  const auto [variant_type, span] = edit.variant->modifyInterval();
  return span.lower() == edit.begin and span.upper() == edit.end;

}

bool appendSequence(kgl::DNA5SequenceLinear& target, const kgl::DNA5SequenceLinear& source) {

  return static_cast<bool>(target.append(source));

}

} // namespace

std::expected<kgl::ConsensusSequence, kgl::ConsensusError>
kgl::buildConsensus(const DNA5SequenceLinear& reference,
                    const OpenRightUnsigned& origin,
                    const OffsetAccounting& accounting,
                    const ConsensusOptions& options) {

  if (reference.length() != origin.size()) {

    return std::unexpected(ConsensusError::ReferenceOriginMismatch);

  }

  if (reference.empty()) {

    return std::unexpected(ConsensusError::EmptyReference);

  }

  ConsensusSequence result;
  result.origin_ = origin;

  // 5' clip from upstream deletes: D = [d0,d1) with d0 < R0 < d1.
  ContigOffset_t w0 = origin.lower();
  for (auto const& edit : accounting.edits()) {

    if (edit.kind == EditKind::Delete and edit.begin < origin.lower() and edit.end > w0) {

      w0 = std::max(w0, edit.end);

    }

  }

  // Retained reference span: empty when an upstream delete consumed the whole window.
  if (w0 <= origin.upper()) {

    result.retained_reference_ = {w0, origin.upper()};

  } else {

    result.retained_reference_ = {origin.upper(), origin.upper()};
    w0 = origin.upper();

  }

  // Select the edits that affect the retained window.
  //   Delete : intersects [w0, R1).
  //   SNP    : contained in [w0, R1).
  //   Insert : member span [i0-1, i0+1) contained in [w0, R1), i.e. w0+1 <= i0 <= R1-1.
  std::vector<Edit> applied;
  for (auto const& edit : accounting.edits()) {

    if (edit.kind == EditKind::Delete) {

      if (edit.end > w0 and edit.begin < origin.upper()) {

        applied.push_back(edit);

      }

    } else if (edit.kind == EditKind::Insert) {

      if (edit.begin >= w0 + 1 and edit.begin < origin.upper()) {

        applied.push_back(edit);

      }

    } else {

      if (edit.begin >= w0 and edit.begin < origin.upper()) {

        applied.push_back(edit);

      }

    }

  }

  // Single forward pass over the retained window.
  DNA5SequenceLinear consensus;
  consensus.reserve(result.retained_reference_.size());

  const std::string_view reference_view = reference.getStringView();
  const auto ref_base = [&](ContigOffset_t contig_offset) -> char {
    return reference_view[contig_offset - origin.lower()];
  };

  const auto emit_gap = [&](ContigOffset_t gap_begin, ContigOffset_t gap_end) -> bool {
    if (gap_begin >= gap_end) {

      return true;

    }
    // Rebase contig coordinates to the reference sequence's local [0, n) interval.
    auto gap_opt = reference.subSequence({gap_begin - origin.lower(), gap_end - origin.lower()});
    if (not gap_opt) {

      return false;

    }
    return appendSequence(consensus, gap_opt.value());

  };

  ContigOffset_t read = w0;
  for (auto const& edit : applied) {

    if (edit.begin < read) {

      continue;   // shadowed by a preceding delete (defensive; should not occur)

    }

    if (not emit_gap(read, std::min(edit.begin, origin.upper()))) {

      return std::unexpected(ConsensusError::LengthInvariant);

    }

    switch (edit.kind) {

      case EditKind::Snp: {
        if (edit.begin >= origin.upper()) {

          break;

        }
        if (options.verify_reference
            and ref_base(edit.begin) != edit.variant->reference().getStringView().front()) {

          return std::unexpected(ConsensusError::ReferenceMismatch);

        }
        consensus.push_back(edit.variant->alternate()[0]);
        read = edit.begin + 1;
      }
      break;

      case EditKind::Insert: {
        if (edit.begin < w0 + 1 or edit.begin >= origin.upper()) {

          break;

        }
        if (options.verify_reference
            and ref_base(edit.begin - 1) != edit.variant->reference().getStringView().front()) {

          return std::unexpected(ConsensusError::ReferenceMismatch);

        }
        for (const char nucleotide : edit.insertedPayload()) {

          consensus.push_back(DNA5::convertChar(nucleotide));

        }
        read = edit.begin;
      }
      break;

      case EditKind::Delete: {
        const ContigOffset_t delete_begin = std::max(edit.begin, w0);
        const ContigOffset_t delete_end = std::min(edit.end, origin.upper());
        if (delete_begin < delete_end and options.verify_reference and deleteVariantMatches(edit)) {

          const ContigSize_t length = delete_end - delete_begin;
          const std::string_view reference_allele = edit.variant->reference().getStringView();
          // Canonical reference is anchor + deleted bases. The index of contig base x
          // within the allele is 1 + (x - edit.begin); for a clipped upstream delete
          // delete_begin may exceed edit.begin.
          const ContigSize_t offset_in_allele =
              edit.variant->referenceSize() - edit.consumed() + (delete_begin - edit.begin);
          const std::string_view expected = reference_allele.substr(offset_in_allele, length);
          auto actual_opt = reference.subSequence({delete_begin - origin.lower(), delete_end - origin.lower()});
          if (not actual_opt or actual_opt.value().getStringView() != expected) {

            return std::unexpected(ConsensusError::ReferenceMismatch);

          }

        }
        read = delete_end;
      }
      break;

    }

  }

  // Tail.
  if (not emit_gap(std::min(read, origin.upper()), origin.upper())) {

    return std::unexpected(ConsensusError::LengthInvariant);

  }

  result.length_delta_ = static_cast<SignedOffset_t>(consensus.length())
                         - static_cast<SignedOffset_t>(result.retained_reference_.size());
  result.modified_interval_ = {w0, w0 + consensus.length()};
  result.bases_ = std::move(consensus);
  result.applied_ = std::move(applied);

  if (result.bases_.length() != result.modified_interval_.size()) {

    return std::unexpected(ConsensusError::LengthInvariant);

  }

  return result;

}

std::expected<kgl::ConsensusSequence, kgl::ConsensusError>
kgl::buildConsensus(const DNA5SequenceLinear& sequence,
                    const OffsetAccounting& accounting,
                    const ConsensusOptions& options) {

  return buildConsensus(sequence, sequence.interval(), accounting, options);

}

std::optional<kgl::DNA5SequenceLinear>
kgl::ConsensusSequence::slice(const OpenRightUnsigned& reference_interval) const {

  if (not origin_.containsInterval(reference_interval)) {

    return std::nullopt;

  }

  // Local index immediately before reference base x, for x in [w0, R1].
  //   inserts : the payload precedes base i0, so x >= i0 includes it
  //   deletes : x inside [b, e] collapses to the deletion boundary
  //   SNPs    : the alternate base occupies the single consumed position
  const auto local_of = [this](ContigOffset_t x) -> ContigOffset_t {

    if (x <= retained_reference_.lower()) {

      return 0;   // 5'-clipped/upstream region collapses to the retained start

    }

    ContigOffset_t local{0};
    ContigOffset_t read = retained_reference_.lower();

    for (auto const& edit : applied_) {

      if (edit.begin > x) {

        break;

      }

      if (edit.begin > read) {

        local += edit.begin - read;
        read = edit.begin;

      }

      switch (edit.kind) {

        case EditKind::Insert:
          // Boundary convention (matches the reference mapping): an offset exactly at
          // the insert point maps to the payload start, so an interval beginning at the
          // point includes the payload. Offsets past the point include it via the shift.
          if (edit.begin == x) {

            return local;

          }
          local += edit.inserted();
          break;

        case EditKind::Snp:
          if (edit.begin == x) {

            return local;

          }
          local += 1;
          read = edit.end;
          break;

        case EditKind::Delete:
          if (x <= edit.end) {

            return local;
          }
          read = edit.end;
          break;

      }

    }

    return local + (x - read);

  };

  const ContigOffset_t lower = local_of(reference_interval.lower());
  const ContigOffset_t upper = local_of(reference_interval.upper());

  if (upper <= lower) {

    return DNA5SequenceLinear{};   // wholly deleted

  }

  return bases_.subSequence({lower, upper});

}

bool kgl::ConsensusSequence::isDeleted(const OpenRightUnsigned& reference_interval) const {

  const auto slice_opt = slice(reference_interval);
  return slice_opt and slice_opt->empty();

}
