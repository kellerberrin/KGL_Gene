//
// Created by kellerberrin on 22/11/20.
//

#ifndef KGL_KGL_ANALYSIS_MUTATION_SYNTHETIC_H
#define KGL_KGL_ANALYSIS_MUTATION_SYNTHETIC_H



#include "kgl_ped_parser.h"
#include "kgl_analysis_mutation_inbreed_calc.h"
#include "kgl_analysis_mutation_inbreed_output.h"


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
  static bool syntheticInbreeding( std::shared_ptr<const PopulationVariant> unphased_ptr,
                                   const InbreedingParameters& paramaters,
                                   InbreedingOutputResults& results);

private:

  // Synthetic inbreeding values.
  constexpr static const double MIN_INBREEDING_COEFICIENT = -0.5; // Set the minimum inbreeding coefficient
  constexpr static const double MAX_INBREEDING_COEFICIENT = 0.5; // Set the maximum inbreeding coefficient
  constexpr static const double STEP_INBREEDING_COEFICIENT = 0.01; // Set the inbreeding step

  [[nodiscard]] static ResultsMap processSynResults( std::shared_ptr<const PopulationVariant> unphased_ptr,
                                                     const InbreedingParameters& parameters);

};



} // namespace




#endif //KGL_KGL_ANALYSIS_MUTATION_SYNTHETIC_H
