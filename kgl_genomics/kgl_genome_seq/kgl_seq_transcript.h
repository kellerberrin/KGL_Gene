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

  SequenceTranscript(const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                     std::shared_ptr<const TranscriptionSequence> transcript_ptr )
                     : transcript_ptr_(std::move(transcript_ptr)) {

    const auto [statistics, status] = createModifiedSequence( contig_variant_ptr);
    sequence_stats_ = statistics;
    sequence_status_ = status;

  }
  ~SequenceTranscript() = default;

  // This can be called multiple times on the same object. An updated AdjustedSequence object is created each time.


  // Returns a sequence of the concatenated and modified exons. Not in strand sense.
  [[nodiscard]] std::optional<DNA5SequenceLinear> getModifiedGene() const;

  // Returns a sequence of the concatenated and original unmodified exons. Not in strand sense.
  [[nodiscard]] std::optional<DNA5SequenceLinear> getOriginalGene() const;

  // The adjusted sequence object has the original interval, detailed internal sequence structure and s
  // modified and original sequences.
  [[nodiscard]] const AdjustedSequence& adjustedSequence() const { return adjusted_sequence_; }
  [[nodiscard]] const SequenceStats& sequenceStatistics() const { return sequence_stats_; }
  [[nodiscard]] bool sequenceStatus() const { return sequence_status_; }

private:

  AdjustedSequence adjusted_sequence_;
  std::shared_ptr<const TranscriptionSequence> transcript_ptr_;
  SequenceStats sequence_stats_;
  bool sequence_status_;

  [[nodiscard]] std::pair<SequenceStats, bool> createModifiedSequence(const std::shared_ptr<const ContigDB>& contig_variant_ptr);

};



} // namespace

#endif //KGL_SEQ_TRANSCRIPT_H
