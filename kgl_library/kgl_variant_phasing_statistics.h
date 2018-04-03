//
// Created by kellerberrin on 3/04/18.
//

#ifndef KGL_VARIANT_PHASING_STATISTICS_H
#define KGL_VARIANT_PHASING_STATISTICS_H

#include "kgl_variant_factory_vcf_phasing.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object generates phasing statistics. In particular, overlapping variant statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


class PopulationPhasingStatistics {

public:

  explicit PopulationPhasingStatistics() = default;
  PopulationPhasingStatistics(const PopulationPhasingStatistics&) = default;
  ~PopulationPhasingStatistics() = default;

  bool phasedSNPs(const VCFPopulation& vcf_population,
                  std::shared_ptr<const GenomeDatabase> genome_db);

private:


};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_PHASING_STATISTICS_H
