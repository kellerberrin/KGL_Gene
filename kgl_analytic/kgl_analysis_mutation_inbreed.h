//
// Created by kellerberrin on 23/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_H
#define KGL_ANALYSIS_MUTATION_INBREED_H


#include "kgl_variant_db_phased.h"
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

  // Analyze population and synthetic data for all defined inbreeding algorithms.
  static bool InbreedingAll( std::shared_ptr<const GenomeReference> genome_ptr,
                             std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                             std::shared_ptr<const DiploidPopulation> diploid_population,
                             std::shared_ptr<const GenomePEDData> ped_data,
                             InbreedingParameters& paramaters,
                             InbreedingOutputResults& results);

private:

  // Analyze both population and synthetic data for a specified algorithm.
  [[nodiscard]] static bool Inbreeding( std::shared_ptr<const GenomeReference> genome_ptr,
                                        std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                        std::shared_ptr<const DiploidPopulation> diploid_population,
                                        std::shared_ptr<const GenomePEDData> ped_data,
                                        InbreedingParameters& paramaters,
                                        InbreedingOutputResults& results);


  // Analyze a presented diploid population for inbreeding.
  [[nodiscard]] static bool populationInbreeding(std::shared_ptr<const GenomeReference> genome_ptr,
                                                 std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                 const DiploidPopulation& diploid_population,
                                                 const GenomePEDData& ped_data,
                                                 InbreedingParameters& paramaters,
                                                 InbreedingOutputResults& results);

  // Analyze a presented diploid population for inbreeding.
  [[nodiscard]] static bool singleInbreeding(std::shared_ptr<const GenomeReference> genome_ptr,
                                                 std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                 const DiploidPopulation& diploid_population,
                                                 const GenomePEDData& ped_data,
                                                 InbreedingParameters& paramaters,
                                                 InbreedingOutputResults& results);

  [[nodiscard]] static ResultsMap populationInbreedingSample( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                              const DiploidPopulation& diploid_population,
                                                              const GenomePEDData& ped_data,
                                                              const InbreedingParameters& parameters);

  // Use a threadpool to calculate the inbreeding coefficients.
  [[nodiscard]] static ResultsMap processResults( const ContigLocusMap& contig_locus_map,
                                                  const DiploidPopulation& diploid_population,
                                                  const GenomePEDData& ped_data,
                                                  const InbreedingParameters& parameters);



};



} // namespace





#endif //KGL_ANALYSIS_MUTATION_INBREED_H
