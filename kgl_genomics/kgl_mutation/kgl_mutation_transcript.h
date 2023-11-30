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

  // Returns a sequence of the concatenated and reference unmodified exons. Not in strand sense.
  [[nodiscard]] std::optional<DNA5SequenceLinear> getOriginalLinear() const;

  // In strand sense. Returns a sequence of the concatenated and modified exons.
  [[nodiscard]] std::optional<DNA5SequenceCoding> getModifiedCoding() const;

  // In strand sense. Returns a sequence of the concatenated and reference unmodified exons.
  [[nodiscard]] std::optional<DNA5SequenceCoding> getOriginalCoding() const;

  // In strand sense. Returns a sequence of the concatenated and modified exons.
  // The coding sequence is also analysed for protein validity (ncRNA just return 'NCRNA').
  // The size_t returns the size of the amino sequence including the first stop sequence.
  [[nodiscard]] std::optional<std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t>> getModifiedValidity() const;
  // Adjust for mod3.
  [[nodiscard]] std::optional<std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t>> getModifiedAdjustedValidity() const;

  // In strand sense. Returns a sequence of the concatenated and reference unmodified exons.
  // The coding sequence is also analysed for protein validity (ncRNA just return 'NCRNA').
  // The size_t returns the size of the amino sequence including the first stop sequence.
  [[nodiscard]] std::optional<std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t>> getOriginalValidity() const;





  // The adjusted sequence object has the reference interval, detailed internal sequence structure and s
  // modified and reference sequences.
  [[nodiscard]] const AdjustedSequence& adjustedSequence() const { return adjusted_sequence_; }
  [[nodiscard]] const FilteredVariantStats& filterStatistics() const { return filter_stats_; }
  [[nodiscard]] bool sequenceStatus() const { return sequence_status_; }

private:

  AdjustedSequence adjusted_sequence_;
  std::shared_ptr<const TranscriptionSequence> transcript_ptr_;
  FilteredVariantStats filter_stats_;
  bool sequence_status_;
  constexpr static const ContigOffset_t PRIME_3_BUFFER_{200} ;
  OpenRightUnsigned prime_3_extend_{0, 0};   // The returned extension may not be the requested extension.
  constexpr static const ContigOffset_t PRIME_5_BUFFER_{0} ;
  OpenRightUnsigned prime_5_extend_{0, 0};   // The returned extension may not be the requested extension.

  [[nodiscard]] std::pair<FilteredVariantStats, bool>
  createModifiedSequence(const std::shared_ptr<const ContigDB>& contig_variant_ptr, SeqVariantFilterType filter_type);
  // Try to adjust the coding sequence to create a viable protein sequence (adjust mod3 and no stop).
  [[nodiscard]] std::optional<DNA5SequenceLinear> getModifiedAdjusted() const;
  [[nodiscard]] std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t> getValidity(DNA5SequenceLinear&& linear_coding) const;



};



} // namespace

#endif //KGL_SEQ_TRANSCRIPT_H
