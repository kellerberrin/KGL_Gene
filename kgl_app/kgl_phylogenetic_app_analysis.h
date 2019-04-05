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
                       std::shared_ptr<const GenomeCollection> genome_collection_ptr,
                       std::shared_ptr<const UnphasedPopulation> unphased_population_ptr,
                       std::shared_ptr<const PhasedPopulation> population_ptr) :  runtime_options_(runtime_options),
                                                                                  genome_collection_ptr_(genome_collection_ptr),
                                                                                  unphased_population_ptr_(unphased_population_ptr),
                                                                                  population_ptr_(population_ptr) {

    analysis_map_[ANALYZE_INTERVAL_] = &PhylogeneticAnalysis::performInterval;
    analysis_map_[ANALYZE_SEQUENCES_] = &PhylogeneticAnalysis::performSequence;
    analysis_map_[ANALYZE_REGION_] = &PhylogeneticAnalysis::performRegion;
    analysis_map_[ANALYZE_UPGMA_] = &PhylogeneticAnalysis::performUPGMA;
    analysis_map_[ANALYZE_GENE_] = &PhylogeneticAnalysis::performGene;
    analysis_map_[ANALYZE_RNA_] = &PhylogeneticAnalysis::performRNA;
    analysis_map_[ANALYZE_SNP_] = &PhylogeneticAnalysis::performSNP;
    analysis_map_[ANALYZE_MIX_] = &PhylogeneticAnalysis::performMix;
    analysis_map_[ANALYZE_MOTIF_] = &PhylogeneticAnalysis::performMotif;
    analysis_map_[ANALYZE_FINESTRUCTURE_] = &PhylogeneticAnalysis::performFineStructure;

  }
  ~PhylogeneticAnalysis() = default;

  // Analysis dispatcher
  void performAnalysis(const std::string& analysis_type);


private:

  // Dynamic analysis dispatcher map
  AnalysisMap analysis_map_;

  // Analysis data and options.
  const RuntimeProperties& runtime_options_;
  std::shared_ptr<const GenomeCollection> genome_collection_ptr_;
  std::shared_ptr<const UnphasedPopulation> unphased_population_ptr_;
  std::shared_ptr<const PhasedPopulation> population_ptr_;

  // Analytic types.
  static constexpr const char* ANALYZE_INTERVAL_ = "INTERVAL";
  static constexpr const char* ANALYZE_SEQUENCES_ = "SEQUENCE";
  static constexpr const char* ANALYZE_REGION_ = "REGION";
  static constexpr const char* ANALYZE_UPGMA_ = "UPGMA";
  static constexpr const char* ANALYZE_GENE_ = "GENE";
  static constexpr const char* ANALYZE_RNA_ = "RNA";
  static constexpr const char* ANALYZE_SNP_ = "SNP";
  static constexpr const char* ANALYZE_MIX_ = "MIX";
  static constexpr const char* ANALYZE_MOTIF_ = "MOTIF";
  static constexpr const char* ANALYZE_FINESTRUCTURE_ = "FINE_STRUCTURE";

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
  void performFineStructure();

  std::shared_ptr<const GenomeDatabase> get3D7Genome() const { return genome_collection_ptr_->get3D7Genome(); }
  std::shared_ptr<const GenomeDatabase> getHB3Genome() const { return genome_collection_ptr_->getHB3Genome(); }
  std::shared_ptr<const GenomeDatabase> getITGenome() const { return genome_collection_ptr_->getITGenome(); }
  std::shared_ptr<const GenomeDatabase> getDD2Genome() const { return genome_collection_ptr_->getDD2Genome(); }
  std::shared_ptr<const GenomeDatabase> get7GBGenome() const { return genome_collection_ptr_->get7GBGenome(); }
  std::shared_ptr<const GenomeDatabase> firstGenome() const { return genome_collection_ptr_->getMap().begin()->second; }

};




}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_PHYLOGENETIC_APP_ANALYSIS_H
