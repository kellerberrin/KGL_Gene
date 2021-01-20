//
// Created by kellerberrin on 7/12/20.
//

#include "kgl_analysis_inbreed_execute.h"
#include "kgl_analysis_inbreed_diploid.h"
#include "kgl_analysis_inbreed_synthetic.h"
#include "kgl_variant_filter.h"


namespace kgl = kellerberrin::genome;



// Perform the genetic analysis per iteration.
bool kgl::ExecuteInbreedingAnalysis::executeAnalysis(std::shared_ptr<const PopulationDB> diploid_population,
                                                     std::shared_ptr<const PopulationDB> unphased_population,
                                                     std::shared_ptr<const GenomePEDData> ped_data,
                                                     InbreedParamOutput& parameters) {


  if (parameters.getParameters().analyzeSynthetic()) {

    // Check that we have what we need.
    if (not unphased_population) {

      ExecEnv::log().error("InbreedAnalysis::iterationAnalysis; Insufficient data, cannot process synthetic diploid inbreeding");
      return false;

    }

    return processSynthetic(unphased_population, parameters);

  } else {

    // Check that we have what we need.
    if (not (diploid_population and unphased_population)) {

      ExecEnv::log().error("ExecuteInbreedingAnalysis::processDiploid; Insufficient data, cannot process diploid inbreeding");
      return false;

    }

    return processDiploid(diploid_population, unphased_population, ped_data, parameters);

  }

}


// Perform the genetic analysis per iteration.
bool kgl::ExecuteInbreedingAnalysis::processDiploid(std::shared_ptr<const PopulationDB> diploid_population,
                                                    std::shared_ptr<const PopulationDB> unphased_population,
                                                    std::shared_ptr<const GenomePEDData> ped_data,
                                                    InbreedParamOutput& parameters) {


  InbreedingAnalysis::populationInbreeding( unphased_population, *diploid_population, *ped_data, parameters);

  return true;

}

// Perform an inbreeding analysis of a synthetic population.
bool kgl::ExecuteInbreedingAnalysis::processSynthetic(std::shared_ptr<const PopulationDB> unphased_population,
                                                      InbreedParamOutput& parameters) {

  SyntheticAnalysis::syntheticInbreeding(unphased_population, parameters);

  return true;

}

