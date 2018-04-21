//
// Created by kellerberrin on 21/04/18.
//

#ifndef KGL_PHYLOGENETIC_APP_ANALYSIS_H
#define KGL_PHYLOGENETIC_APP_ANALYSIS_H


#include "kgl_genome_db.h"
#include "kgl_vcf_parser_data.h"
#include "kgl_phylogenetic_env.h"

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
                              std::shared_ptr<const ParserAnalysis> parser_analysis_ptr);


};




}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_PHYLOGENETIC_APP_ANALYSIS_H
