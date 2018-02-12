//
// Created by kellerberrin on 12/02/18.
//

#ifndef KGL_PHLOGENETIC_APP_H
#define KGL_PHLOGENETIC_APP_H

#include "kgl_phylogenetic_env.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



// Application class implements the mainline logic and controls
// data object lifetimes, see kgl_phylogenetic_app.cc.
class PhylogeneticApp {

public:

  PhylogeneticApp(const Phylogenetic& args);
  ~PhylogeneticApp() = default;

private:

  std::shared_ptr<const GenomeVariant>  getGenomeVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                          const std::string& file_name,
                                                          bool vcf_is_gatk,
                                                          const std::string& genome_name,
                                                          Phred_t read_quality,
                                                          Phred_t variant_quality,
                                                          long min_count,
                                                          double min_proportion);

  void performAnalysis(const Phylogenetic& args,
                       std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                       std::shared_ptr<const PopulationVariant> pop_variant_ptr);

};




} //  organization level namespace
}  // project level namespace


#endif // KGL_PHLOGENETIC_APP_H
