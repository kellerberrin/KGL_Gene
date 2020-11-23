//
// Created by kellerberrin on 24/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_AUX_H
#define KGL_ANALYSIS_MUTATION_INBREED_AUX_H

#include "kgl_variant_db_phased.h"
#include "kgl_analysis_mutation_inbreed_freqdb.h"

#include <memory>
#include <map>

namespace kellerberrin::genome {   //  organization::project level namespace


// Allele locus lists indexed by superpopulation
using LocusMap = std::map<std::string, std::shared_ptr<const ContigVariant>>;
// LocusMaps indexed by contig id.
using ContigLocusMap = std::map<ContigId_t, LocusMap>;

// Just a namespace.
class InbreedSampling  {

public:

  InbreedSampling() = delete;
  ~InbreedSampling() = delete;


  // Uses the defined contigs in the unphased population to create a contig map of population locuii.
  [[nodiscard]] static ContigLocusMap getPopulationLocusMap(  std::shared_ptr<const UnphasedPopulation> population_ptr,
                                                              double min_af,
                                                              double max_af,
                                                              ContigOffset_t locii_spacing,
                                                              ContigOffset_t upper_offset = DEFAULT_UPPER_OFFSET,
                                                              ContigOffset_t lower_offset = DEFAULT_LOWER_OFFSET);

  constexpr static const ContigOffset_t DEFAULT_UPPER_OFFSET = 1000000000; // Larger than all contigs
  constexpr static const ContigOffset_t DEFAULT_LOWER_OFFSET = 0;

private:


  // Get a list of potential allele locus with a specified spacing to minimise linkage dis-equilibrium
  // and at a specified frequency for the super population. Used as a template for calculating
  // the inbreeding coefficient and sample relatedness
  using LocusReturnPair = std::pair<std::string, std::shared_ptr<const ContigVariant>>;
  [[nodiscard]] static LocusReturnPair getLocusList( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                     const ContigId_t& contig_id,
                                                     const std::string& super_population,
                                                     double min_af,
                                                     double max_af,
                                                     ContigOffset_t locii_spacing,
                                                     ContigOffset_t upper_offset,
                                                     ContigOffset_t lower_offset);

  // Generate a list of locii to sample a population for the inbreeding coefficient.
  [[nodiscard]] static LocusMap getPopulationLocus(std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                   const ContigId_t& contig_id,
                                                   double min_af,
                                                   double max_af,
                                                   ContigOffset_t locii_spacing,
                                                   ContigOffset_t upper_offset,
                                                   ContigOffset_t lower_offset);

};





} // namespace

#endif //KGL_KGL_ANALYSIS_MUTATION_INBREED_AUX_H
