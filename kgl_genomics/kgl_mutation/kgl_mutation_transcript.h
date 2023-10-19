//
// Created by kellerberrin on 16/09/23.
//

#ifndef KGL_SEQ_TRANSCRIPT_H
#define KGL_SEQ_TRANSCRIPT_H

#include "kgl_mutation_sequence.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class SequenceTranscript {

public:


  SequenceTranscript(const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                     std::shared_ptr<const TranscriptionSequence> transcript_ptr,
                     SeqVariantFilterType filter_type = SeqVariantFilterType::DEFAULT_SEQ_FILTER)
                     : transcript_ptr_(std::move(transcript_ptr)) {

    const auto [statistics, status] = createModifiedSequence( contig_variant_ptr, filter_type);
    filter_stats_ = statistics;
    sequence_status_ = status;

  }
  ~SequenceTranscript() = default;

  // This can be called multiple times on the same object. An updated AdjustedSequence object is created each time.


  // Returns a sequence of the concatenated and modified exons. Not in strand sense.
  [[nodiscard]] std::optional<DNA5SequenceLinear> getModifiedLinear() const;

  // Returns a sequence of the concatenated and original unmodified exons. Not in strand sense.
  [[nodiscard]] std::optional<DNA5SequenceLinear> getOriginalLinear() const;

  // In strand sense. Returns a sequence of the concatenated and modified exons.
  [[nodiscard]] std::optional<DNA5SequenceCoding> getModifiedCoding() const;

  // In strand sense. Returns a sequence of the concatenated and original unmodified exons.
  [[nodiscard]] std::optional<DNA5SequenceCoding> getOriginalCoding() const;

  // In strand sense. Returns a sequence of the concatenated and modified exons.
  // The coding sequence is also analysed for protein validity (ncRNA just return 'NCRNA').
  // The size_t returns the size of the amino sequence including the first stop sequence.
  [[nodiscard]] std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t> getModifiedValidity() const;

  // In strand sense. Returns a sequence of the concatenated and original unmodified exons.
  // The coding sequence is also analysed for protein validity (ncRNA just return 'NCRNA').
  // The size_t returns the size of the amino sequence including the first stop sequence.
  [[nodiscard]] std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t> getOriginalValidity() const;


  // The adjusted sequence object has the original interval, detailed internal sequence structure and s
  // modified and original sequences.
  [[nodiscard]] const AdjustedSequence& adjustedSequence() const { return adjusted_sequence_; }
  [[nodiscard]] const FilteredVariantStats& filterStatistics() const { return filter_stats_; }
  [[nodiscard]] bool sequenceStatus() const { return sequence_status_; }

private:

  AdjustedSequence adjusted_sequence_;
  std::shared_ptr<const TranscriptionSequence> transcript_ptr_;
  FilteredVariantStats filter_stats_;
  bool sequence_status_;

  [[nodiscard]] std::pair<FilteredVariantStats, bool>
  createModifiedSequence(const std::shared_ptr<const ContigDB>& contig_variant_ptr, SeqVariantFilterType filter_type);


};



} // namespace

#endif //KGL_SEQ_TRANSCRIPT_H
