//
// Created by kellerberrin on 21/04/18.
//

#ifndef KGL_PHYLOGENETIC_APP_ANALYSIS_H
#define KGL_PHYLOGENETIC_APP_ANALYSIS_H


#include "kgl_genome_db.h"
#include "kgl_variant_db.h"
#include "kgl_phylogenetic_app.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level analysis object.
// Dispatches analytics.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PhylogeneticAnalysis {

public:

using AnalysisFuncPtr = void (PhylogeneticAnalysis::*)();
using AnalysisMap = std::map<std::string, AnalysisFuncPtr>;


  PhylogeneticAnalysis(const RuntimeProperties& runtime_options,
                       std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                       std::shared_ptr<const UnphasedPopulation> unphased_population_ptr,
                       std::shared_ptr<const PhasedPopulation> population_ptr) :  runtime_options_(runtime_options),
                                                                                  genome_db_ptr_(genome_db_ptr),
                                                                                  unphased_population_ptr_(unphased_population_ptr),
                                                                                  population_ptr_(population_ptr) {

    analysis_map_[ANALYZE_INTERVAL] = &PhylogeneticAnalysis::performInterval;
    analysis_map_[ANALYZE_SEQUENCES] = &PhylogeneticAnalysis::performSequence;
    analysis_map_[ANALYZE_REGION] = &PhylogeneticAnalysis::performRegion;
    analysis_map_[ANALYZE_UPGMA] = &PhylogeneticAnalysis::performUPGMA;
    analysis_map_[ANALYZE_GENE] = &PhylogeneticAnalysis::performGene;
    analysis_map_[ANALYZE_RNA] = &PhylogeneticAnalysis::performRNA;
    analysis_map_[ANALYZE_SNP] = &PhylogeneticAnalysis::performSNP;
    analysis_map_[ANALYZE_MIX] = &PhylogeneticAnalysis::performMix;
    analysis_map_[ANALYZE_MOTIF] = &PhylogeneticAnalysis::performMotif;

  }
  ~PhylogeneticAnalysis() = default;

  // Analysis dispatcher
  void performAnalysis(const std::string& analysis_type);


private:

  // Dynamic analysis dispatcher map
  AnalysisMap analysis_map_;

  // Analysis data and options.
  const RuntimeProperties& runtime_options_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  std::shared_ptr<const UnphasedPopulation> unphased_population_ptr_;
  std::shared_ptr<const PhasedPopulation> population_ptr_;

  // Analytic types.
  static constexpr const char* ANALYZE_INTERVAL = "INTERVAL";
  static constexpr const char* ANALYZE_SEQUENCES = "SEQUENCE";
  static constexpr const char* ANALYZE_REGION = "REGION";
  static constexpr const char* ANALYZE_UPGMA = "UPGMA";
  static constexpr const char* ANALYZE_GENE = "GENE";
  static constexpr const char* ANALYZE_RNA = "RNA";
  static constexpr const char* ANALYZE_SNP = "SNP";
  static constexpr const char* ANALYZE_MIX = "MIX";
  static constexpr const char* ANALYZE_MOTIF = "MOTIF";

  // Analysis routines
  void performMotif();
  void performSequence();
  void performInterval();
  void performGene();
  void performRegion();
  void performMix();
  void performSNP();
  void performUPGMA();
  void performRNA();

};




}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_PHYLOGENETIC_APP_ANALYSIS_H
