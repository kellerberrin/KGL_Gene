//
// Created by kellerberrin on 17/09/23.
//

#ifndef KGL_MUTATION_TRANSLATE_H
#define KGL_MUTATION_TRANSLATE_H

#include "kel_interval.h"
#include "kgl_genome_types.h"
#include "kgl_mutation_variant_map.h"

#include <map>


namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These helper classes keep track of the offset translation between the zero-offset modified sequence
// and the contig based offset.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AdjustedModifiedOffset {

public:

  AdjustedModifiedOffset(ContigOffset_t contig_offset,
                         ContigOffset_t modified_offset,
                         SignedOffset_t offset_adjust,
                         const SequenceVariantUpdate& sequence_update)
      : contig_offset_(contig_offset),
        modified_offset_(modified_offset),
        offset_adjust_(offset_adjust),
        sequence_update_(sequence_update) {}

  ~AdjustedModifiedOffset() = default;

  [[nodiscard]] ContigOffset_t contigOffset() const { return contig_offset_; }

  [[nodiscard]] ContigOffset_t modifiedOffset() const { return modified_offset_; }

  [[nodiscard]] SignedOffset_t cumulativeAdjust() const { return cumulative_adjust_; }

  [[nodiscard]] SignedOffset_t offsetAdjust() const { return offset_adjust_; }

  [[nodiscard]] const SequenceVariantUpdate& sequenceUpdate() const { return sequence_update_; }


  [[nodiscard]] std::string toString() const;

  // This is the updated cumulative indel size/offset adjustment used in the offset map below.
  void updateCumulative(SignedOffset_t cumulative) { cumulative_adjust_ = cumulative; }


private:

  ContigOffset_t contig_offset_;
  ContigOffset_t modified_offset_;
  SignedOffset_t offset_adjust_;
  SignedOffset_t cumulative_adjust_{0};
  SequenceVariantUpdate sequence_update_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using AdjustOffsetMap = std::map<ContigOffset_t, AdjustedModifiedOffset>;

class ModifiedOffsetMap {

public:

  explicit ModifiedOffsetMap(const OpenRightUnsigned &contig_interval) : contig_interval_(contig_interval) {}
  ~ModifiedOffsetMap() = default;

  // Add an indel offset adjust to the map.
  // Note that the actual indel increment is passed as an argument.
  // But, importantly, the cumulative indel offset is calculated and also stored in the map.
  [[nodiscard]] bool addModifiedOffset(AdjustedModifiedOffset modified_offset);

  // Given a contig based sequence interval, return the equivalent zero-based indel modified interval (adjusts for delete indels).
  [[nodiscard]] std::pair<OpenRightUnsigned, bool> lookupModifiedInterval(const OpenRightUnsigned &contig_interval) const;

  // Given a contig based sequence interval, return the equivalent zero-based unmodified/original interval.
  [[nodiscard]] std::pair<OpenRightUnsigned, bool> lookupOriginalInterval(const OpenRightUnsigned &contig_interval) const;

  // Given an interval map offset, returns the equivalent offset into the zero-offset modified sequence [0, n).
  // Calculates the offset into the modified_sequence_.
  [[nodiscard]] std::pair<ContigOffset_t, bool> modifiedZeroOffset(ContigOffset_t offset) const;

  // Given an interval map offset, returns the equivalent offset into the zero-offset original sequence [0, m), where (n-m) = base_offset_adjust_.
  // Calculates the offset into the original_sequence_.
  [[nodiscard]] std::pair<ContigOffset_t, bool> originalZeroOffset(ContigOffset_t contig_offset) const { return calcModifiedOffset(contig_offset, 0); }

  [[nodiscard]] const OpenRightUnsigned& contigInterval() const { return contig_interval_; }

  void reinitialize(const OpenRightUnsigned &contig_interval) { contig_interval_ = contig_interval; adjust_offset_map_.clear(); }

private:

  OpenRightUnsigned contig_interval_;
  AdjustOffsetMap adjust_offset_map_;

  // Lookup an indel modified zero-offset sequence - generally a gene or similar.
  // Note that we must account for contig offsets that occur in the shadow of a delete.
  [[nodiscard]] std::pair<ContigOffset_t, bool> lookupIndelOffset(ContigOffset_t offset) const;

  // Get the previous indel adjustment to the supplied contig offset.
  [[nodiscard]] std::optional<AdjustedModifiedOffset> getPreviousIndel(ContigOffset_t contig_offset) const;
  // Given a contig offset and an indel adjustment, calculate a zero-based offset - this does not adjust for delete indels.
  [[nodiscard]] std::pair<ContigOffset_t, bool> calcModifiedOffset(ContigOffset_t contig_offset, SignedOffset_t indel_adjust) const;
  // Finds the map offset element that is < contig_offset. If no such element exists then returns 0.
  [[nodiscard]] SignedOffset_t getIndelOffset(ContigOffset_t contig_offset) const;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////




} // namespace



#endif //KGL_MUTATION_TRANSLATE_H
