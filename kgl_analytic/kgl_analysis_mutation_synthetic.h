//
// Created by kellerberrin on 22/11/20.
//

#ifndef KGL_KGL_ANALYSIS_MUTATION_SYNTHETIC_H
#define KGL_KGL_ANALYSIS_MUTATION_SYNTHETIC_H



#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"
#include "kgl_analysis_mutation_inbreed_calc.h"


namespace kellerberrin::genome {   //  organization::project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Performs the inbreeding analysis.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class SyntheticAnalysis  {

public:

  SyntheticAnalysis() = delete;
  ~SyntheticAnalysis() = delete;

  // Construct a synthetic population and analyze it.
  // The synthetic population is constructed from the unphased population.
  [[nodiscard]] static bool syntheticInbreeding( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                 InbreedingParameters& paramaters);

private:

  // Synthetic inbreeding values.
  constexpr static const double MIN_INBREEDING_COEFICIENT = -0.5; // Set the minimum inbreeding coefficient
  constexpr static const double MAX_INBREEDING_COEFICIENT = 0.5; // Set the maximum inbreeding coefficient
  constexpr static const double STEP_INBREEDING_COEFICIENT = 0.01; // Set the inbreeding step



  bool static syntheticInbreedingSample( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                         InbreedingParameters& parameters);

  [[nodiscard]] static ResultsMap processSynResults( const ContigLocusMap& contig_locus_map,
                                                     InbreedingParameters& parameters);

};



} // namespace




#endif //KGL_KGL_ANALYSIS_MUTATION_SYNTHETIC_H
