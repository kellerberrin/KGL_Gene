//
// Created by kellerberrin on 7/09/23.
//

#ifndef KGL_MUTATION_SEQUENCE_H
#define KGL_MUTATION_SEQUENCE_H


#include "kel_interval_unsigned.h"
#include "kgl_mutation_translate.h"
#include "kgl_mutation_variant_map.h"
#include "kgl_genome_contig.h"


namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This class actually modifies the zero-based sequence.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Alias to reduce code noise
using DualSeqOpt = std::optional<std::pair<DNA5SequenceLinear, DNA5SequenceLinear>>;

class AdjustedSequence {

public:

  AdjustedSequence() :  modified_offset_map_({0, 0}) { clear(); }
  ~AdjustedSequence() = default;

  // The coding interval for the contig_ref_ptr sequence.
  [[nodiscard]] const OpenRightUnsigned& contigInterval() const { return modified_offset_map_.contigInterval(); }
  // The modify map for the zero based modified sequence.
  [[nodiscard]] const ModifiedOffsetMap& sequenceMap() const { return modified_offset_map_; }
  // Flag set if we have valid sequences.
  [[nodiscard]] bool validModifiedSequence() const { return  valid_modified_sequence_; }
  // The zero-based modified sub-sequence.
  // Note that the specified sub_interval must be contained in the initial contig_interval
  [[nodiscard]] std::optional<DNA5SequenceLinear> modifiedSubSequence(const OpenRightUnsigned& sub_interval) const;
  // The zero-based unmodified sub-sequence.
  // Note that the specified sub_interval must be contained in the initial contig_interval
  [[nodiscard]] std::optional<DNA5SequenceLinear> originalSubSequence(const OpenRightUnsigned& sub_interval) const;

  // Update the unmodified zero-based sequence into the modified sequence.
  // This function can update the modified sequence multiple times.
  [[nodiscard]] bool updateSequence(const std::shared_ptr<const ContigReference>& contig_ref_ptr,
                                    const SequenceVariantFilter& filtered_variants);

  // This function moves the original (.first) and modified (.second) sequences and initializes (clears) the object.
  DualSeqOpt moveSequenceClear();

private:

  IntervalModifyMap interval_modify_map_;
  ModifiedOffsetMap modified_offset_map_;
  DNA5SequenceLinear modified_sequence_;
  DNA5SequenceLinear original_sequence_;
  bool valid_modified_sequence_{false};

  // Update the unmodified zero-based sequence into the modified sequence.
  // This function can update the modified sequence multiple times.
  [[nodiscard]] bool updateSequence(const std::shared_ptr<const ContigReference>& contig_ref_ptr,
                                    const OpenRightUnsigned& contig_interval,
                                    const IntervalModifyMap& interval_modify_map);

  void initializeSequences(const std::shared_ptr<const ContigReference>& contig_ref_ptr,
                           const OpenRightUnsigned& contig_interval,
                           const IntervalModifyMap& interval_modify_map);

  // Update the unmodified zero-based sequence into the modified sequence.
  [[nodiscard]] bool updateSequence();
  // Reset all data structures including sequences, validModifiedSequence() set to false.
  void clear();
  // Verify updated sequence.
  [[nodiscard]] bool verifyUpdatedSequence() const;

  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceSNP(const SequenceVariantUpdate& sequence_update);
  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceDelete(const SequenceVariantUpdate& sequence_update);
  [[nodiscard]] std::pair<bool, SignedOffset_t> updateSequenceInsert(const SequenceVariantUpdate& sequence_update);

};



}  // Namespace



#endif //KGL_MUTATION_SEQUENCE_H
