//
// Created by kellerberrin on 22/11/20.
//

#ifndef KGL_KGA_ANALYSIS_INBREED_SYNTHETIC_H
#define KGL_KGA_ANALYSIS_INBREED_SYNTHETIC_H



#include "kgl_hsgenealogy_parser.h"
#include "kga_analysis_inbreed_calc.h"
#include "kga_analysis_inbreed_output.h"


namespace kellerberrin::genome::analysis {   //  organization::project level namespace



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
  static bool syntheticInbreeding(std::shared_ptr<const PopulationDB> unphased_ptr, InbreedParamOutput& param_output);

private:

  // Synthetic inbreeding values.
  constexpr static const double MIN_INBREEDING_COEFICIENT = -0.5; // Set the minimum inbreeding coefficient
  constexpr static const double MAX_INBREEDING_COEFICIENT = 0.5; // Set the maximum inbreeding coefficient
  constexpr static const double STEP_INBREEDING_COEFICIENT = 0.01; // Set the inbreeding step

  [[nodiscard]] static ResultsMap processSynResults( std::shared_ptr<const PopulationDB> unphased_ptr,
                                                     const InbreedingParameters& parameters);

};



} // namespace




#endif //KGL_KGA_ANALYSIS_INBREED_SYNTHETIC_H
