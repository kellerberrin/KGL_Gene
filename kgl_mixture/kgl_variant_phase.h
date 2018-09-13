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

  static bool fileHaploidPhasing(const std::string& phase_file,
                                 size_t vcf_ploidy,
                                 std::shared_ptr<const UnphasedPopulation> vcf_population_ptr,
                                 std::shared_ptr<const GenomeDatabase> genome_db,
                                 std::shared_ptr<PhasedPopulation> haploid_population);

  // A simple phasing strategy that removes all conflicting variants.
  static bool haploidPhasing(size_t vcf_ploidy,
                             std::shared_ptr<const UnphasedPopulation> vcf_population_ptr,
                             std::shared_ptr<const GenomeDatabase> genome_db,
                             std::shared_ptr<PhasedPopulation> haploid_population);

private:

  // Returns false if no heterozygous variant selected, else the selected variant return as an index to the vector.
  static bool analyseCountStatistics(const UnphasedVectorVariantCount& unphased_vector, size_t& phase_index);

  // The proportion required for a heterozygous variant to be accepted as the haploid variant.
  constexpr static const double HETEROZYGOUS_PROPORTION_ = 0.9;

};






}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_PHASE_H
