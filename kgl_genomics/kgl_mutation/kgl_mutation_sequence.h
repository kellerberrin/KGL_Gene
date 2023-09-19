//
// Created by kellerberrin on 7/09/23.
//

#ifndef KGL_MUTATION_SEQUENCE_H
#define KGL_MUTATION_SEQUENCE_H


#include "kel_interval.h"
#include "kgl_mutation_translate.h"
#include "kgl_mutation_variant_map.h"
#include "kgl_genome_contig.h"


namespace kellerberrin::genome {   //  organization::project level namespace

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
  // The zero-based modified sub-sequence.
  // Note that the specified sub_interval must be contained in the initial contig_interval
  [[nodiscard]] std::optional<DNA5SequenceLinear> modifiedSubSequence(const OpenRightUnsigned& sub_interval) const;
  // The zero-based unmodified sub-sequence.
  // Note that the specified sub_interval must be contained in the initial contig_interval
  [[nodiscard]] std::optional<DNA5SequenceLinear> originalSubSequence(const OpenRightUnsigned& sub_interval) const;
  // Update the unmodified zero-based sequence into the modified sequence.
  [[nodiscard]] bool updateSequence();

private:

  std::shared_ptr<const ContigReference> contig_ref_ptr_;
  const OpenRightUnsigned contig_interval_;
  IntervalModifyMap interval_modify_map_;
  ModifiedOffsetMap modified_offset_map_;
  DNA5SequenceLinear modified_sequence_;
  DNA5SequenceLinear original_sequence_;

  void initializeSequences();

  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceSNP(const SequenceVariantUpdate& sequence_update);
  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceDelete(const SequenceVariantUpdate& sequence_update);
  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceInsert(const SequenceVariantUpdate& sequence_update);

};



}  // Namespace



#endif //KGL_MUTATION_SEQUENCE_H
