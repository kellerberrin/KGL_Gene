//
// Created by kellerberrin on 23/04/18.
//

#ifndef KGL_VARIANT_PHASE_H
#define KGL_VARIANT_PHASE_H



#include "kel_utility.h"
#include "kgl_variant_db.h"


namespace kellerberrin::genome {   //  organization::project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object accepts unphased variants from the VCF parser and phases them.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


class GenomePhasing {

public:

  explicit GenomePhasing() = default;
  ~GenomePhasing() = default;

  static std::shared_ptr<UnphasedPopulation> filterClonal(const std::string& phase_file,
                                                          std::shared_ptr<const UnphasedPopulation> unphased_population_ptr);

  // A simple phasing strategy that removes all conflicting variants.
  static bool haploidPhasing(size_t vcf_ploidy,
                             const std::shared_ptr<const UnphasedPopulation>& vcf_population_ptr,
                             const std::shared_ptr<HaploidPopulation>& haploid_population);

private:

  // Returns false if no heterozygous variant selected, else the selected variant return as an index to the vector.
  static bool analyseCountStatistics(const VirtualContigOffset& unphased_vector, size_t& phase_index);

  // The proportion required for a heterozygous variant to be accepted as the haploid variant.
  constexpr static const double HETEROZYGOUS_PROPORTION_ = 0.9;

};






}   // end namespace


#endif //KGL_VARIANT_PHASE_H
