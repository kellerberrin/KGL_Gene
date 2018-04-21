//
// Created by kellerberrin on 21/04/18.
//

#ifndef KGL_PHYLOGENETIC_APP_ANALYSIS_H
#define KGL_PHYLOGENETIC_APP_ANALYSIS_H


#include "kgl_genome_db.h"
#include "kgl_phylogenetic_env.h"
#include "kgl_variant_factory_vcf_phasing.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



// The top level analysis object.
class PhylogeneticAnalysis {

public:

  PhylogeneticAnalysis() = delete; // Singleton
  ~PhylogeneticAnalysis() = delete;

  static bool checkAnalysisType(const std::string& AnalysisType);

  static void performAnalysis(const Phylogenetic& args,
                              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                              std::shared_ptr<const VCFPopulation> vcf_population_ptr);
private:

  // Analytic types.
  static constexpr const char* ANALYZE_INTERVAL = "INTERVAL";
  static constexpr const char* ANALYZE_SEQUENCES = "SEQUENCE";
  static constexpr const char* ANALYZE_REGION = "REGION";
  static constexpr const char* ANALYZE_UPGMA = "UPGMA";
  static constexpr const char* ANALYZE_GENE = "GENE";
  static constexpr const char* ANALYZE_RNA = "RNA";
  static constexpr const char* ANALYZE_SNP = "SNP";

};




}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_PHYLOGENETIC_APP_ANALYSIS_H
