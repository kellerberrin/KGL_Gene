//
// Created by kellerberrin on 21/04/18.
//

#ifndef KGL_PHYLOGENETIC_APP_ANALYSIS_H
#define KGL_PHYLOGENETIC_APP_ANALYSIS_H


#include "kgl_variant_db_population.h"
#include "kgl_gene_app.h"

namespace kellerberrin::genome {   //  organization::project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Legacy P. falciparum code. All functionality to be shifted to packages.
// Leave for the moment to ensure the code still compiles.
// The top level analysis object.
// Dispatches analytics.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PhylogeneticAnalysis {

public:

using AnalysisFuncPtr = void (PhylogeneticAnalysis::*)();
using AnalysisMap = std::map<std::string, AnalysisFuncPtr>;


  PhylogeneticAnalysis(const RuntimeProperties& runtime_options,
                       std::shared_ptr<const GenomeCollection> genome_collection_ptr,
                       std::shared_ptr<const PopulationDB> unphased_population_ptr,
                       std::shared_ptr<const PopulationDB> population_ptr) : runtime_options_(runtime_options),
                                                                             genome_collection_ptr_(genome_collection_ptr),
                                                                             unphased_population_ptr_(unphased_population_ptr),
                                                                             population_ptr_(population_ptr) {

    analysis_map_[ANALYZE_SEQUENCES_] = &PhylogeneticAnalysis::performSequence;
    analysis_map_[ANALYZE_GENE_] = &PhylogeneticAnalysis::performGene;
    analysis_map_[ANALYZE_RNA_] = &PhylogeneticAnalysis::performRNA;
    analysis_map_[ANALYZE_SNP_] = &PhylogeneticAnalysis::performSNP;
    analysis_map_[ANALYZE_MIX_] = &PhylogeneticAnalysis::performMix;
    analysis_map_[ANALYZE_MOTIF_] = &PhylogeneticAnalysis::performMotif;

  }
  ~PhylogeneticAnalysis() = default;

  // Simple analysis parser. Perform a list a analysis types
  // Delimited by one of the following characters ",;-".
  // For example "GENE;INTERVAL;UPGMA" performs the 3 analysis types.
  void performAnalysis(const std::string& analysis_type);

private:

  static constexpr char ANALYSIS_SEPARATORS_ = ',';
  // Analysis dispatcher. Perform each individual analysis
  void dispatchAnalysis(const std::string& analysis_type);

  // Dynamic analysis dispatcher map
  AnalysisMap analysis_map_;

  // Analysis data and options.
  const RuntimeProperties& runtime_options_;
  std::shared_ptr<const GenomeCollection> genome_collection_ptr_;
  std::shared_ptr<const PopulationDB> unphased_population_ptr_;
  std::shared_ptr<const PopulationDB> population_ptr_;
  std::shared_ptr<const PopulationDB> reference_genome_ptr_;

  // Analytic types.
  static constexpr const char* ANALYZE_SEQUENCES_ = "SEQUENCE";
  static constexpr const char* ANALYZE_REGION_ = "REGION";
  static constexpr const char* ANALYZE_GENE_ = "GENE";
  static constexpr const char* ANALYZE_RNA_ = "RNA";
  static constexpr const char* ANALYZE_SNP_ = "SNP";
  static constexpr const char* ANALYZE_MIX_ = "MIX";
  static constexpr const char* ANALYZE_MOTIF_ = "MOTIF";

  // Analysis routines that are called by text.
  void performMotif();  // MOTIF
  void performSequence(); // SEQUENCE
  void performGene(); // GENE
  void performMix(); // MIX
  void performSNP(); // SNP
  void performRNA(); // RNA



};


}   // end namespace


#endif //KGL_PHYLOGENETIC_APP_ANALYSIS_H
