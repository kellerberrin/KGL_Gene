//
// Created by kellerberrin on 7/09/23.
//

#ifndef KGL_MUTATION_SEQUENCE_H
#define KGL_MUTATION_SEQUENCE_H

#include "kgl_mutation_variant_map.h"
#include "kgl_genome_contig.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class AdjustedSequence {

public:

  AdjustedSequence(std::shared_ptr<const ContigReference> contig_ref_ptr,
                   OpenRightUnsigned contig_interval,
                   IntervalModifyMap interval_modify_map)
                    : contig_ref_ptr_(std::move(contig_ref_ptr)),
                      contig_interval_(contig_interval),
                      base_offset_adjust_(contig_interval_.lower()),
                      interval_modify_map_(std::move(interval_modify_map)) { initializeSequences(); }
  ~AdjustedSequence() = default;

  [[nodiscard]] OpenRightUnsigned contigInterval() const { return contig_interval_; }
  [[nodiscard]] const IntervalModifyMap& intervalMap() const { return interval_modify_map_; }
  [[nodiscard]] const DNA5SequenceLinear& modifiedSequence() const { return modified_sequence_; }
  [[nodiscard]] const DNA5SequenceLinear& originalSequence() const { return original_sequence_; }

  [[nodiscard]] bool updateSequence();

private:


  std::shared_ptr<const ContigReference> contig_ref_ptr_;
  const OpenRightUnsigned contig_interval_;
  const ContigOffset_t base_offset_adjust_; // Initialized to contig_interval_.lower();
  IntervalModifyMap interval_modify_map_;
  IntervalModifyMultiMap update_audit_map_;
  DNA5SequenceLinear modified_sequence_;
  DNA5SequenceLinear original_sequence_;

  // Used to adjust variant and modify map offsets to offsets within modified_sequence_.
  ContigOffset_t current_offset_adjust_{0}; // Updated to priorInterval.lower() foe each update record.
  SignedOffset_t modify_offset_adjust_{0}; // Adjusted by indels, +ve for insert, -ve for delete.

  void initializeSequences();

  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceSNP(const SequenceVariantUpdate& sequence_update);
  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceDelete(const SequenceVariantUpdate& sequence_update);
  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceInsert(const SequenceVariantUpdate& sequence_update);

  // Given an interval map offset, returns the equivalent offset into the zero-offset modified sequence [0, n).
  [[nodiscard]] ContigOffset_t translateZeroOffset(ContigOffset_t offset) const; // Calculates the offset into the modified_sequence_.
  [[nodiscard]] OpenRightUnsigned translateZeroOffset(const OpenRightUnsigned& interval) const; // Same as above but translates an interval.

  // Uses a base that is updated for each modify record.
  ContigOffset_t translateRelativeOffset(ContigOffset_t offset) const;
  OpenRightUnsigned translateRelativeOffset(const OpenRightUnsigned& interval) const;

  // Given an interval map offset, returns the equivalent offset into the zero-offset original sequence [0, m), where (n-m) = base_offset_adjust_.
  [[nodiscard]] ContigOffset_t originalZeroOffset(ContigOffset_t offset) const;  // Calculates the offset into the original_sequence_.
  [[nodiscard]] OpenRightUnsigned translateOriginalOffset(const OpenRightUnsigned& interval) const; // Same as above but translates an interval.

};



}  // Namespace



#endif //KGL_MUTATION_SEQUENCE_H
