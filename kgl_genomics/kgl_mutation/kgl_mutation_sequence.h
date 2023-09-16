//
// Created by kellerberrin on 7/09/23.
//

#ifndef KGL_MUTATION_SEQUENCE_H
#define KGL_MUTATION_SEQUENCE_H

#include "kgl_mutation_variant_map.h"
#include "kgl_genome_contig.h"


namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These helper classes keep track of the offset between the zero-offset modified sequence and the contig based offset.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AdjustedModifiedOffset {

public:

  AdjustedModifiedOffset(ContigOffset_t contig_offset,
                         ContigOffset_t modified_offset,
                         SignedOffset_t indel_adjust)
                         : contig_offset_(contig_offset),
                           modified_offset_(modified_offset),
                           indel_adjust_(indel_adjust) {}
  ~AdjustedModifiedOffset() = default;

  [[nodiscard]] ContigOffset_t contigOffset() const { return contig_offset_; }
  [[nodiscard]] ContigOffset_t modifiedOffset() const { return modified_offset_; }
  [[nodiscard]] SignedOffset_t indelAdjust() const { return indel_adjust_; }

  // This is the updated cumulative indel size/offset adjustment used in the offset map below.
  void addCumulativeOffset(SignedOffset_t update_adjust) { indel_adjust_ += update_adjust; }

private:

  ContigOffset_t contig_offset_;
  ContigOffset_t modified_offset_;
  SignedOffset_t indel_adjust_;

};

using AdjustOffsetMap = std::map<ContigOffset_t, AdjustedModifiedOffset>;

class ModifiedOffsetMap {

public:

  ModifiedOffsetMap(const OpenRightUnsigned& contig_interval) : contig_interval_(contig_interval) {}
  ~ModifiedOffsetMap() = default;

  // Add an indel offset adjust to the map.
  // Note that the actual indel increment is passed as an argument.
  // But, importantly, the cumulative indel offset is stored in the map.
  [[nodiscard]] bool addModifiedOffset(AdjustedModifiedOffset modified_offset);
  // Given a contig based sequence interval, return the equivalent zero-based indel modified interval.
  [[nodiscard]] OpenRightUnsigned convertContigModified(const OpenRightUnsigned& contig_interval) const;
  // Finds the last cumulative offset or 0 if empty map.
  [[nodiscard]] SignedOffset_t getCumulativeOffset() const;

  // Given an interval map offset, returns the equivalent offset into the zero-offset modified sequence [0, n).
  // Calculates the offset into the modified_sequence_.
  [[nodiscard]] std::pair<ContigOffset_t, bool> translateZeroOffset(ContigOffset_t offset) const;

private:

  const OpenRightUnsigned contig_interval_;
  AdjustOffsetMap adjust_offset_map_;

  // Finds the map offset element that is <= contig_offset. If not such element exists then returns 0.
  [[nodiscard]] SignedOffset_t getIndelOffset(ContigOffset_t contig_offset) const;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This class actually modifies the zero-based sequence.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AdjustedSequence {

public:

  AdjustedSequence(std::shared_ptr<const ContigReference> contig_ref_ptr,
                   OpenRightUnsigned contig_interval,
                   IntervalModifyMap interval_modify_map)
                    : contig_ref_ptr_(std::move(contig_ref_ptr)),
                      contig_interval_(contig_interval),
                      interval_modify_map_(std::move(interval_modify_map)),
                      modified_offset_map_(contig_interval) { initializeSequences(); }
  ~AdjustedSequence() = default;

  // The coding interval for the contig sequence.
  [[nodiscard]] const OpenRightUnsigned& contigInterval() const { return contig_interval_; }
  // The modify map for the original unmodified sequence.
  [[nodiscard]] const IntervalModifyMap& intervalMap() const { return interval_modify_map_; }
  // The modify map for the zero based modified sequence.
  [[nodiscard]] const ModifiedOffsetMap& sequenceMap() const { return modified_offset_map_; }
  // The zero-based modified sequence.
  [[nodiscard]] const DNA5SequenceLinear& modifiedSequence() const { return modified_sequence_; }
  // The zero-based unmodified sequence.
  [[nodiscard]] const DNA5SequenceLinear& originalSequence() const { return original_sequence_; }
  // Update the unmodified zero-based sequence into the modified sequence.
  [[nodiscard]] bool updateSequence();

private:

  std::shared_ptr<const ContigReference> contig_ref_ptr_;
  const OpenRightUnsigned contig_interval_;
  IntervalModifyMap interval_modify_map_;
  ModifiedOffsetMap modified_offset_map_;
  DNA5SequenceLinear modified_sequence_;
  DNA5SequenceLinear original_sequence_;

  // Used to adjust variant and modify map offsets to offsets within modified_sequence_.
  SignedOffset_t modify_offset_adjust_{0}; // Adjusted by indels, +ve for insert, -ve for delete.

  void initializeSequences();

  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceSNP(const SequenceVariantUpdate& sequence_update);
  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceDelete(const SequenceVariantUpdate& sequence_update);
  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceInsert(const SequenceVariantUpdate& sequence_update);

  // Given an interval map offset, returns the equivalent offset into the zero-offset modified sequence [0, n).
  [[nodiscard]] std::pair<ContigOffset_t, bool> translateZeroOffset(ContigOffset_t offset) const; // Calculates the offset into the modified_sequence_.

  // Given an interval map offset, returns the equivalent offset into the zero-offset original sequence [0, m), where (n-m) = base_offset_adjust_.
  [[nodiscard]] std::pair<ContigOffset_t, bool> originalZeroOffset(ContigOffset_t offset) const;  // Calculates the offset into the original_sequence_.

};



}  // Namespace



#endif //KGL_MUTATION_SEQUENCE_H
