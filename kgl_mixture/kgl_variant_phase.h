//
// Created by kellerberrin on 23/04/18.
//

#ifndef KGL_VARIANT_PHASE_H
#define KGL_VARIANT_PHASE_H



#include "kgl_utility.h"
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object accepts unphased variants from the VCF parser and phases them.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


class GenomePhasing {

public:

  explicit GenomePhasing() = default;
  ~GenomePhasing() = default;

  static bool haploidPhasing(std::shared_ptr<const UnphasedPopulation> vcf_population_ptr,
                             std::shared_ptr<const GenomeDatabase> genome_db,
                             std::shared_ptr<PhasedPopulation> haploid_population);

private:



};






}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_PHASE_H
