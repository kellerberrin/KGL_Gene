//
// Created by kellerberrin on 17/11/23.
//

#ifndef KGA_ANALYSIS_LIB_SEQ_GENE_H
#define KGA_ANALYSIS_LIB_SEQ_GENE_H

#include "kgl_variant_db_population.h"
#include "kgl_genome_genome.h"

namespace kellerberrin::genome::analysis {   //  organization::project level namespace

struct TranscriptSequenceRecord {
  std::string sequence_string_;
  std::set<GenomeId_t> genomes_;

};


class AnalysisTranscriptSequence {

public:

  AnalysisTranscriptSequence(std::shared_ptr<const TranscriptionSequence>, std::string ident_work_directory);
  ~AnalysisTranscriptSequence() = default;

  void performTranscriptAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr);

  void printReport() const;

private:


  std::shared_ptr<const TranscriptionSequence> transcript_ptr_;
  std::map<std::string, TranscriptSequenceRecord> sequence_map_;

};



} // Namespace

#endif //KGA_ANLYSIS_LIB_SEQ_GENE_H
