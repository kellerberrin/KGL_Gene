//
// Created by kellerberrin on 7/12/20.
//

#ifndef KGL_ANALYSIS_INBREED_EXECUTE_H
#define KGL_ANALYSIS_INBREED_EXECUTE_H


#include "kgl_hsgenealogy_parser.h"
#include "kga_analysis_inbreed_calc.h"
#include "kga_analysis_inbreed_args.h"
#include "kga_analysis_inbreed_output.h"


namespace kellerberrin::genome {   //  organization::project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Sets up the Inbreeding Parameters and then executes the analysis.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ExecuteInbreedingAnalysis {

public:

  // The population variant data.
  ExecuteInbreedingAnalysis() = delete;
  ~ExecuteInbreedingAnalysis() = delete;

  static bool executeAnalysis(std::shared_ptr<const PopulationDB> diploid_population,
                              std::shared_ptr<const PopulationDB> unphased_population,
                              std::shared_ptr<const HsGenomeGenealogyData> ped_data,
                              InbreedParamOutput& param_output);

private:

  static bool processDiploid(std::shared_ptr<const PopulationDB> diploid_population,
                             std::shared_ptr<const PopulationDB> unphased_population,
                             std::shared_ptr<const HsGenomeGenealogyData> ped_data,
                             InbreedParamOutput& param_output);
  static bool processSynthetic(std::shared_ptr<const PopulationDB> unphased_population,
                               InbreedParamOutput& param_output);

};


} // namespace




#endif //KGL_ANALYSIS_INBREED_EXECUTE_H
