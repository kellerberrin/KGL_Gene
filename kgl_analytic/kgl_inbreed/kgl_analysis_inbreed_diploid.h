//
// Created by kellerberrin on 23/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_H
#define KGL_ANALYSIS_MUTATION_INBREED_H


#include "kgl_genealogy_parser.h"
#include "kgl_analysis_inbreed_calc.h"
#include "kgl_analysis_inbreed_output.h"


namespace kellerberrin::genome {   //  organization::project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Performs the inbreeding analysis.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class InbreedingAnalysis  {

public:

  InbreedingAnalysis() = delete;
  ~InbreedingAnalysis() = delete;


  // Analyze a presented diploid population for inbreeding.
  static bool populationInbreeding(std::shared_ptr<const PopulationDB> unphased_ptr,
                                   const PopulationDB& diploid_population,
                                   const GenomeGenealogyData& ped_data,
                                   InbreedParamOutput& param_output);

private:

  [[nodiscard]] static ResultsMap populationInbreedingSample( std::shared_ptr<const PopulationDB> unphased_ptr,
                                                              const PopulationDB& diploid_population,
                                                              const GenomeGenealogyData& ped_data,
                                                              const InbreedingParameters& parameters);

  // Use a threadpool to calculate the inbreeding coefficients.
  [[nodiscard]] static ResultsMap processResults( const ContigLocusMap& contig_locus_map,
                                                  const PopulationDB& diploid_population,
                                                  const GenomeGenealogyData& ped_data,
                                                  const InbreedingParameters& parameters);

};



} // namespace





#endif //KGL_ANALYSIS_MUTATION_INBREED_H
