//
// Created by kellerberrin on 16/09/23.
//

#ifndef KGL_SEQ_TRANSCRIPT_H
#define KGL_SEQ_TRANSCRIPT_H

#include "kgl_mutation_sequence.h"
#include "kgl_seq_interval.h"


namespace kellerberrin::genome {   //  organization::project level namespace

// Simple struct to return generated sequence stats.
struct SequenceStats {

  size_t map_size_{0};
  size_t non_unique_count_{0};
  size_t upstream_deleted_{0};
  CodingSequenceValidity original_sequence_{CodingSequenceValidity::VALID_PROTEIN};
  CodingSequenceValidity modified_sequence_{CodingSequenceValidity::VALID_PROTEIN};

};


class SequenceTranscript {

public:

  SequenceTranscript() = default;
  ~SequenceTranscript() = default;

  // This can be called multiple times on the same object. An updated AdjustedSequence object is created each time.
  [[nodiscard]] std::pair<SequenceStats, bool> createModifiedSequence( const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                                                                       const std::shared_ptr<const ContigReference>& contig_reference_ptr,
                                                                       const OpenRightUnsigned& modified_interval);


  // Returns a sequence of the concatenated and modified exons. Not in strand sense.
  [[nodiscard]] std::optional<DNA5SequenceLinear> getModifiedGene( const GeneIntervalStructure& gene_interval,
                                                                   const FeatureIdent_t& transcript_id) const;

  // Returns a sequence of the concatenated and original unmodified exons. Not in strand sense.
  [[nodiscard]] std::optional<DNA5SequenceLinear> getOriginalGene( const GeneIntervalStructure& gene_interval,
                                                                   const FeatureIdent_t& transcript_id) const;

  // The adjusted sequence object has the original interval, detailed internal sequence structure and s
  // modified and original sequences.
  [[nodiscard]] const AdjustedSequence& adjustedSequence() const { return adjusted_sequence_; }

private:

  AdjustedSequence adjusted_sequence_;

  // Given a vector of intervals, sorted them by lower() and then retrieves the modified sequence,
  // and concatenates these in sorted sequence order.
  [[nodiscard]] std::optional<DNA5SequenceLinear> concatModifiedSequences(const std::vector<OpenRightUnsigned>& interval_vector) const ;

  // Given a vector of intervals, sorted them by lower() and then retrieves the original sequence,
  // and concatenates these in sorted sequence order.
  [[nodiscard]] std::optional<DNA5SequenceLinear> concatOriginalSequences(const std::vector<OpenRightUnsigned>& interval_vector) const;


};



} // namespace

#endif //KGL_SEQ_TRANSCRIPT_H
