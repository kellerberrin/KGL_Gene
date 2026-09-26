//
// Consensus library (glm-refactor polish of deepseek-refactor).
//
// Module 4: gene/transcript client.
// Keeps the exact public surface of the reference SequenceTranscript so that all
// client translation units compile unchanged. All exon, strand, codon and validity
// logic lives here; modules 1-3 know nothing about genes.
//

#ifndef KGL_CONSENSUS_TRANSCRIPT_H
#define KGL_CONSENSUS_TRANSCRIPT_H

#include "kgl_consensus_select.h"
#include "kgl_consensus_sequence.h"
#include "kgl_genome_prelim.h"

#include <optional>
#include <tuple>

namespace kellerberrin::genome {

// Alias retained for the clients' includes and for the compatibility facade.
using DualSeqOpt = std::optional<std::pair<DNA5SequenceLinear, DNA5SequenceLinear>>;

class SequenceTranscript {

public:

  SequenceTranscript(const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                     std::shared_ptr<const TranscriptionSequence> transcript_ptr,
                     SeqVariantFilterType filter_type = SeqVariantFilterType::DEFAULT_SEQ_FILTER);
  ~SequenceTranscript() = default;

  SequenceTranscript(const SequenceTranscript&) = delete;
  SequenceTranscript& operator=(const SequenceTranscript&) = delete;

  // Returns a sequence of the concatenated and modified exons. Not in strand sense.
  [[nodiscard]] std::optional<DNA5SequenceLinear> getModifiedLinear() const;

  // Returns a sequence of the concatenated and reference unmodified exons. Not in strand sense.
  [[nodiscard]] std::optional<DNA5SequenceLinear> getOriginalLinear() const;

  // In strand sense. Returns a sequence of the concatenated and modified exons.
  [[nodiscard]] std::optional<DNA5SequenceCoding> getModifiedCoding() const;

  // In strand sense. Returns a sequence of the concatenated and reference unmodified exons.
  [[nodiscard]] std::optional<DNA5SequenceCoding> getOriginalCoding() const;

  // In strand sense. The coding sequence is analysed for protein validity.
  // The size_t returns the amino sequence size including the first stop codon.
  [[nodiscard]] std::optional<std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t>> getModifiedValidity() const;
  [[nodiscard]] std::optional<std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t>> getModifiedAdjustedValidity() const;
  [[nodiscard]] std::optional<std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t>> getOriginalValidity() const;

  [[nodiscard]] const FilteredVariantStats& filterStatistics() const { return filter_stats_; }
  [[nodiscard]] bool sequenceStatus() const { return sequence_status_; }

private:

  std::shared_ptr<const TranscriptionSequence> transcript_ptr_;
  FilteredVariantStats filter_stats_;
  bool sequence_status_{false};

  // The consensus over the extended transcript window (module 3 output).
  std::optional<ConsensusSequence> consensus_;
  // Original reference over the same window, for originalSubSequence().
  DNA5SequenceLinear original_window_;
  OpenRightUnsigned window_{0, 0};

  constexpr static const ContigOffset_t PRIME_3_BUFFER_{200};
  constexpr static const ContigOffset_t PRIME_5_BUFFER_{0};
  OpenRightUnsigned prime_3_extend_{0, 0};
  OpenRightUnsigned prime_5_extend_{0, 0};

  void createModifiedSequence(const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                              SeqVariantFilterType filter_type);

  [[nodiscard]] std::optional<DNA5SequenceLinear> modifiedLinear() const;
  [[nodiscard]] std::optional<DNA5SequenceLinear> originalLinear() const;
  [[nodiscard]] std::optional<DNA5SequenceLinear> getModifiedAdjusted() const;
  [[nodiscard]] std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t>
  getValidity(DNA5SequenceLinear&& linear_coding) const;

};

} // namespace

#endif //KGL_CONSENSUS_TRANSCRIPT_H
