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



// The top level analysis object.
class PhylogeneticAnalysis {

public:

  PhylogeneticAnalysis() = delete; // Singleton
  ~PhylogeneticAnalysis() = delete;

  static bool checkAnalysisType(const std::string& AnalysisType);

  static void performAnalysis(const Phylogenetic& args,
                              const RuntimeOptions& runtime_options,
                              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                              std::shared_ptr<const UnphasedPopulation> unphased_population_ptr,
                              std::shared_ptr<const PhasedPopulation> population_ptr);
private:

  // Analytic types.
  static constexpr const char* ANALYZE_INTERVAL = "INTERVAL";
  static constexpr const char* ANALYZE_SEQUENCES = "SEQUENCE";
  static constexpr const char* ANALYZE_REGION = "REGION";
  static constexpr const char* ANALYZE_UPGMA = "UPGMA";
  static constexpr const char* ANALYZE_GENE = "GENE";
  static constexpr const char* ANALYZE_RNA = "RNA";
  static constexpr const char* ANALYZE_SNP = "SNP";
  static constexpr const char* ANALYZE_MIX = "MIX";


};




}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_PHYLOGENETIC_APP_ANALYSIS_H
