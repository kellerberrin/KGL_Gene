//
// Created by kellerberrin on 24/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_AUX_H
#define KGL_ANALYSIS_MUTATION_INBREED_AUX_H

#include "kgl_analysis_inbreed_args.h"
#include "kgl_analysis_inbreed_freq.h"

#include <memory>
#include <map>

namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object uses a range of selection criteria to generate a vector of locii offsets used to select locii for
// inbreeding sampling.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class RetrieveLociiVector {

public:

  RetrieveLociiVector() =delete;
  ~RetrieveLociiVector() = delete;

  static std::vector<ContigOffset_t> getLociiFromTo(std::shared_ptr<const ContigDB> unphased_contig_ptr,
                                                    const std::string& super_population,
                                                    const LociiVectorArguments& arguments);

  static std::vector<ContigOffset_t> getLociiCount(std::shared_ptr<const ContigDB> unphased_contig_ptr,
                                                   const std::string& super_population,
                                                   const LociiVectorArguments& arguments);

private:


  static std::vector<AlleleFreqVector> getAllelesFromTo( std::shared_ptr<const ContigDB> unphased_contig_ptr,
                                                         const std::string& super_population,
                                                         const LociiVectorArguments& arguments);

  static std::vector<AlleleFreqVector> getAllelesCount( std::shared_ptr<const ContigDB> unphased_contig_ptr,
                                                        const std::string& super_population,
                                                        const LociiVectorArguments& arguments);

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object actually uses the locii vector generated above to package the selected locii into a ContigLocusMap
// for presentation to the inbreeding analytics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Allele locus lists indexed by superpopulation
using LocusMap = std::map<std::string, std::shared_ptr<const ContigDB>>;
// LocusMaps indexed by contig id.
using ContigLocusMap = std::map<ContigId_t, LocusMap>;

// Just a namespace.
class InbreedSampling  {

public:

  InbreedSampling() = delete;
  ~InbreedSampling() = delete;


  // Uses the defined contigs in the unphased population to create a contig map of population locii.
  [[nodiscard]] static ContigLocusMap getPopulationLocusMap(  std::shared_ptr<const PopulationDB> population_ptr,
                                                              const LociiVectorArguments& locii_args);

private:


  // Get a list of potential allele locus with a specified spacing to minimise linkage dis-equilibrium
  // and at a specified frequency for the super population. Used as a template for calculating
  // the inbreeding coefficient and sample relatedness
  using LocusReturnPair = std::pair<std::string, std::shared_ptr<const ContigDB>>;
  [[nodiscard]] static LocusReturnPair getLocusList( std::shared_ptr<const PopulationDB> unphased_ptr,
                                                     const ContigId_t& contig_id,
                                                     const std::string& super_population,
                                                     const LociiVectorArguments& locii_args);

  // Generate a list of locii to sample a population for the inbreeding coefficient.
  [[nodiscard]] static LocusMap getPopulationLocus(std::shared_ptr<const PopulationDB> unphased_ptr,
                                                   const ContigId_t& contig_id,
                                                   const LociiVectorArguments& locii_args);

};





} // namespace

#endif //KGL_KGL_ANALYSIS_MUTATION_INBREED_AUX_H
