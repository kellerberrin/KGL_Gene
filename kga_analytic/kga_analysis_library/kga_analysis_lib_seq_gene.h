//
// Created by kellerberrin on 17/11/23.
//

#ifndef KGA_ANALYSIS_LIB_SEQ_GENE_H
#define KGA_ANALYSIS_LIB_SEQ_GENE_H

#include "kgl_variant_db_population.h"
#include "kgl_genome_genome.h"
#include "kgl_mutation_variant_filter.h"

#include "kel_utility.h"

namespace kellerberrin::genome::analysis {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeTuple = std::tuple<DNA5SequenceCoding, CodingSequenceValidity, size_t>;

class TranscriptModifyRecord {

public:

  TranscriptModifyRecord( GenomeId_t genome,
                          std::shared_ptr<const TranscriptionSequence> transcript_ptr,
                          GenomeTuple&& reference,
                          GenomeTuple&& modified)
  : genome_(std::move(genome)),
    transcript_ptr_(std::move(transcript_ptr)),
    reference_transcript_(std::move(reference)),
    modified_transcript_(std::move(modified)) {}
  TranscriptModifyRecord(TranscriptModifyRecord&& rval_copy) noexcept
  : genome_(std::move(rval_copy.genome_)),
    transcript_ptr_(std::move(rval_copy.transcript_ptr_)),
    reference_transcript_(std::move(rval_copy.reference_transcript_)),
    modified_transcript_(std::move(rval_copy.modified_transcript_)) {}
  ~TranscriptModifyRecord() = default;

  [[nodiscard]] const GenomeId_t& genome() const { return genome_; }
  [[nodiscard]] const std::shared_ptr<const TranscriptionSequence>& transcript() const { return transcript_ptr_; }
  [[nodiscard]] const GenomeTuple& reference() const { return reference_transcript_; }
  [[nodiscard]] const GenomeTuple& modified() const { return modified_transcript_; }

private:

  GenomeId_t genome_;
  std::shared_ptr<const TranscriptionSequence> transcript_ptr_;
  GenomeTuple reference_transcript_;
  GenomeTuple modified_transcript_;

};
using GenomeRecordOpt = std::optional<TranscriptModifyRecord>;

struct TranscriptSequenceRecord {

  std::set<GenomeId_t> genomes_;
  std::set<FeatureIdent_t> transcripts_;
  double distance_{0.0};

};
using TranscriptMap = std::map<std::string, TranscriptSequenceRecord>;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  constexpr static const SeqVariantFilterType FILTERTYPE_{SeqVariantFilterType::DEFAULT_SEQ_FILTER};

  [[nodiscard]] GenomeRecordOpt genomeTranscriptMutation(const std::shared_ptr<const GenomeDB>& genome_db_ptr);

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AnalysisTranscriptFamily {

public:

  AnalysisTranscriptFamily() = default;
  ~AnalysisTranscriptFamily() = default;

  void createAnalysisVector(const GeneVector& gene_vector);
  void performFamilyAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr);
  void printAllReports(const std::string& analysis_directory, const std::string& analysis_sub_directory) const;


private:

  std::vector<AnalysisTranscriptSequence> analysis_vector_;

  AnalysisTranscriptSequence generateTotal() const;

};



} // Namespace

#endif //KGA_ANLYSIS_LIB_SEQ_GENE_H
