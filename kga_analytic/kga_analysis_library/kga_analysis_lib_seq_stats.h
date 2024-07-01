//
// Created by kellerberrin on 9/12/23.
//

#include "kgl_distance_tree_upgma.h"
#include "kga_analysis_lib_seq_gene.h"


#ifndef KGA_ANALYSIS_LIB_SEQ_STATS_H
#define KGA_ANALYSIS_LIB_SEQ_STATS_H


namespace kellerberrin::genome::analysis {   //  organization::project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Top level analysis and report object.
using AnalysisTranscriptMap = std::map<FeatureIdent_t, AnalysisTranscriptSequence>;

class AnalysisTranscriptFamily {

public:

  AnalysisTranscriptFamily() = default;
  ~AnalysisTranscriptFamily() = default;

  void createAnalysisVector(const GeneVector& gene_vector);
  void performFamilyAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr);
  void printAllReports(const std::string& analysis_directory,
                       const std::string& analysis_sub_directory,
                       const std::map<GenomeId_t, std::string>& annotation_map) const;

  void createClassificationTree(const std::string& analysis_directory,
                                const std::string& analysis_sub_directory,
                                const std::map<GenomeId_t, std::string>& annotations,
                                const std::function<bool(const std::string&)>& sample_selection) const;

  [[nodiscard]] const AnalysisTranscriptMap& getMap() const { return analysis_map_; }

private:

  AnalysisTranscriptMap analysis_map_;
  [[nodiscard]] AnalysisTranscriptSequence generateTotal() const;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Indexed by genome and then by transcript id.
using FeatureTranscriptMap = std::map<FeatureIdent_t , std::shared_ptr<const TranscriptModifyRecord>>;
using GenomeTranscriptMap = std::map<GenomeId_t, FeatureTranscriptMap>;

class GenomeTranscriptAnalysis {

public:

  GenomeTranscriptAnalysis() = default;
  ~GenomeTranscriptAnalysis() = default;

  [[nodiscard]] bool createTranscriptMap(const AnalysisTranscriptMap& transcript_map);
  void printGenomeReport(const std::string& report_directory, const std::map<GenomeId_t, std::string>& annotations = {}) const;
  bool createClassificationTree(const std::string& newick_directory,
                                const std::string& transcript_id,
                                const std::map<GenomeId_t, std::string>& annotations,
                                const std::function<bool(const std::string&)>& sample_selection);


private:

  GenomeTranscriptMap genome_map_;

  [[nodiscard]] bool processTranscriptMap(const AnalysisTranscriptMap& transcript_map);
  [[nodiscard]] std::optional<std::vector<FeatureIdent_t>> checkGenomeMap() const;

  constexpr static const std::string REPORT_EXT_{".csv"};
  constexpr static const std::string REPORT_FIELD_{","};
  constexpr static const std::string REPORT_PREFIX_{"genome"};

  constexpr static const std::string NEWICK_SUFFIX_{"newick"};
  constexpr static const std::string NEWICK_EXT_{".txt"};

};





}


#endif //KGA_ANALYSIS_LIB_SEQ_STATS_H
