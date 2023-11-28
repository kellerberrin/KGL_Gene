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
using GenomeRecordOpt = std::optional<GenomeTuple>;

struct TranscriptSequenceRecord {

  std::set<GenomeId_t> genomes_;

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
  ~AnalysisTranscriptSequence() = default;

  void performTranscriptAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr);
  void printReport(const std::string& report_directory) const;

private:


  std::shared_ptr<const TranscriptionSequence> transcript_ptr_;
  const std::string report_directory_;
  std::map<std::string, TranscriptSequenceRecord> sequence_map_;

  TranscriptMap transcript_map_;

  constexpr static const std::string REPORT_EXT_{".csv"};
  constexpr static const std::string REPORT_FIELD_{","};
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

};



} // Namespace

#endif //KGA_ANLYSIS_LIB_SEQ_GENE_H
