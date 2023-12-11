//
// Created by kellerberrin on 17/11/23.
//

#ifndef KGA_ANALYSIS_LIB_SEQ_GENE_H
#define KGA_ANALYSIS_LIB_SEQ_GENE_H

#include "kgl_variant_db_population.h"
#include "kgl_genome_genome.h"
#include "kgl_mutation_variant_filter.h"

#include "kel_utility.h"

#include <compare>

namespace kellerberrin::genome::analysis {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class TranscriptModifyRecord {

public:

  TranscriptModifyRecord( GenomeId_t genome,
                          std::shared_ptr<const TranscriptionSequence> transcript_ptr,
                          DNA5SequenceCoding&& reference_coding,
                          CodingSequenceValidity reference_validity,
                          DNA5SequenceCoding&& modified_coding,
                          CodingSequenceValidity modified_validity)
  : genome_(std::move(genome)),
    transcript_ptr_(std::move(transcript_ptr)),
    original_sequence_(std::move(reference_coding)),
    original_validity_(reference_validity),
    modified_sequence_(std::move(modified_coding)),
    modified_validity_(modified_validity) {}
  TranscriptModifyRecord(TranscriptModifyRecord&& rval_copy) noexcept
  : genome_(std::move(rval_copy.genome_)),
    transcript_ptr_(std::move(rval_copy.transcript_ptr_)),
    original_sequence_(std::move(rval_copy.original_sequence_)),
    original_validity_(rval_copy.original_validity_),
    modified_sequence_(std::move(rval_copy.modified_sequence_)),
    modified_validity_(rval_copy.modified_validity_) {}
  ~TranscriptModifyRecord() = default;

  [[nodiscard]] const GenomeId_t& genome() const { return genome_; }
  [[nodiscard]] const std::shared_ptr<const TranscriptionSequence>& transcript() const { return transcript_ptr_; }
  [[nodiscard]] const DNA5SequenceCoding& reference() const { return original_sequence_; }
  [[nodiscard]] CodingSequenceValidity referenceValidity() const { return original_validity_; }
  [[nodiscard]] const DNA5SequenceCoding& modified() const { return modified_sequence_; }
  [[nodiscard]] CodingSequenceValidity modifiedValidity() const { return modified_validity_; }

  // Define an ordering using the spaceship operator. Records are indexed by genome.
  [[nodiscard]] auto operator<=>(const TranscriptModifyRecord &rhs) const { return genome() <=> rhs.genome(); }
  [[nodiscard]] bool operator==(const TranscriptModifyRecord &rhs) const { return genome() == rhs.genome(); }

private:

  GenomeId_t genome_;
  std::shared_ptr<const TranscriptionSequence> transcript_ptr_;
  DNA5SequenceCoding original_sequence_;
  CodingSequenceValidity original_validity_;
  DNA5SequenceCoding modified_sequence_;
  CodingSequenceValidity modified_validity_;

};
using GenomeRecordOpt = std::optional<std::shared_ptr<const TranscriptModifyRecord>>;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using TransSeqGenomeMap = std::multimap<GenomeId_t, std::shared_ptr<const TranscriptModifyRecord>>;
using TransSeqVector = std::vector<std::shared_ptr<const TranscriptModifyRecord>>;
using TransSeqTranscriptMap = std::map<FeatureIdent_t, TransSeqVector>;

class TranscriptSequenceRecord {

public:

  explicit TranscriptSequenceRecord(const std::shared_ptr<const TranscriptModifyRecord>& genome_modify_ptr);
  ~TranscriptSequenceRecord() = default;

  [[nodiscard]] bool addSequenceRecord(const std::shared_ptr<const TranscriptModifyRecord>& genome_modify_ptr);

  [[nodiscard]] const DNA5SequenceCodingView& getModifiedView() const { return modified_sequence_; }
  [[nodiscard]] const TransSeqGenomeMap& getGenomes() const { return genomes_; }
  [[nodiscard]] const TransSeqTranscriptMap& getTranscripts() const { return transcripts_; }

  // Generates a unique label for different sequences. Used as a map index rather than the actual sequence.
  [[nodiscard]] static std::string generateSequenceLabel(const DNA5SequenceCodingView& seq_view);

private:

  const DNA5SequenceCodingView modified_sequence_;
  TransSeqGenomeMap genomes_;
  TransSeqTranscriptMap transcripts_;

  constexpr static const std::string SEQUENCE_LABEL_FORMAT_{"SEQID_{:x}"};

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using TranscriptMap = std::map<std::string, TranscriptSequenceRecord>;

class AnalysisTranscriptSequence {

public:

  explicit AnalysisTranscriptSequence(std::shared_ptr<const TranscriptionSequence> transcript_ptr)
                                     : transcript_ptr_(std::move(transcript_ptr)) {}
  AnalysisTranscriptSequence(const AnalysisTranscriptSequence& copy) = default;
  ~AnalysisTranscriptSequence() = default;

  void performTranscriptAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr);

  void mergeMap(const TranscriptMap& transcript_map);

  void printReport(const std::string& report_directory) const;
  void printReport(const std::string& report_label, const std::string& report_directory) const;

  [[nodiscard]] const TranscriptionSequence& transcript() const { return *transcript_ptr_; }
  [[nodiscard]] const TranscriptMap& transcriptMap() const { return transcript_map_; }

private:

  std::shared_ptr<const TranscriptionSequence> transcript_ptr_;
  TranscriptMap transcript_map_;

  constexpr static const std::string REPORT_EXT_{".csv"};
  constexpr static const std::string REPORT_FIELD_{","};
  constexpr static const std::string REPORT_SUBFIELD_{"-"};
  constexpr static const std::string REPORT_PREFIX_{"transcript"};
  constexpr static const OpenRightUnsigned SUB_VIEW_INTERVAL_{0, 200};
  constexpr static const SeqVariantFilterType FILTERTYPE_{SeqVariantFilterType::DEFAULT_SEQ_FILTER};

  [[nodiscard]] GenomeRecordOpt genomeTranscriptMutation(const std::shared_ptr<const GenomeDB>& genome_db_ptr);

};



} // Namespace

#endif //KGA_ANALYSIS_LIB_SEQ_GENE_H
