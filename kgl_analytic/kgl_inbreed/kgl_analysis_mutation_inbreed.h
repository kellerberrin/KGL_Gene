//
// Created by kellerberrin on 23/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_H
#define KGL_ANALYSIS_MUTATION_INBREED_H


#include "kgl_ped_parser.h"
#include "kgl_analysis_mutation_inbreed_calc.h"
#include "kgl_analysis_mutation_inbreed_output.h"


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
  static bool populationInbreeding( std::shared_ptr<const PopulationVariant> unphased_ptr,
                                    const PopulationVariant& diploid_population,
                                    const GenomePEDData& ped_data,
                                    const InbreedingParameters& paramaters,
                                    InbreedingOutputResults& results);

private:

  [[nodiscard]] static ResultsMap populationInbreedingSample( std::shared_ptr<const PopulationVariant> unphased_ptr,
                                                              const PopulationVariant& diploid_population,
                                                              const GenomePEDData& ped_data,
                                                              const InbreedingParameters& parameters);

  // Use a threadpool to calculate the inbreeding coefficients.
  [[nodiscard]] static ResultsMap processResults( const ContigLocusMap& contig_locus_map,
                                                  const PopulationVariant& diploid_population,
                                                  const GenomePEDData& ped_data,
                                                  const InbreedingParameters& parameters);

};



} // namespace





#endif //KGL_ANALYSIS_MUTATION_INBREED_H
