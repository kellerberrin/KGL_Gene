//
// Created by kellerberrin on 23/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_H
#define KGL_ANALYSIS_MUTATION_INBREED_H


#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"
#include "kgl_analysis_mutation_inbreed_calc.h"


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
                             const std::string& output_file_name,
                             InbreedingParameters& paramaters);

private:

  // Synthetic inbreeding values.
  constexpr static const double MIN_INBREEDING_COEFICIENT = -0.5; // Set the minimum inbreeding coefficient
  constexpr static const double MAX_INBREEDING_COEFICIENT = 0.5; // Set the maximum inbreeding coefficient
  constexpr static const double STEP_INBREEDING_COEFICIENT = 0.01; // Set the inbreeding step


  constexpr static const char DELIMITER_ = ',';
  constexpr static const char* FILE_EXT_ = ".csv";

  using ResultsMap = std::map<GenomeId_t, LocusResults>;

  // Analyze both population and synthetic data for a specified algorithm.
  [[nodiscard]] static bool Inbreeding( std::shared_ptr<const GenomeReference> genome_ptr,
                                        std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                        std::shared_ptr<const DiploidPopulation> diploid_population,
                                        std::shared_ptr<const GenomePEDData> ped_data,
                                        InbreedingParameters& paramaters);

  // Construct a synthetic population and analyze it.
  // The synthetic population is constructed from the unphased population.
  [[nodiscard]] static bool syntheticInbreeding( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                 InbreedingParameters& paramaters);

  // Analyze a presented diploid population for inbreeding.
  [[nodiscard]] static bool populationInbreeding(std::shared_ptr<const GenomeReference> genome_ptr,
                                                 std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                 const DiploidPopulation& diploid_population,
                                                 const GenomePEDData& ped_data,
                                                 InbreedingParameters& paramaters);

  bool static syntheticInbreedingSample( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                         InbreedingParameters& parameters);

  static bool populationInbreedingSample( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                          const DiploidPopulation& diploid_population,
                                          const GenomePEDData& ped_data,
                                          const std::string& output_file_name,
                                          const InbreedingParameters& parameters);

  // Use a threadpool to calculate the inbreeding coefficients.
  static bool processResults( const ContigLocusMap& contig_locus_map,
                              const DiploidPopulation& diploid_population,
                              const GenomePEDData& ped_data,
                              std::ostream& outfile,
                              const InbreedingParameters& parameters);

  static bool processSynResults( const ContigLocusMap& contig_locus_map,
                                 std::ostream& outfile,
                                 InbreedingParameters& parameters);

  // Write the analysis results to a CSV file.
  static bool writeResults( const ContigId_t& contig_id,
                            const ResultsMap& locus_results,
                            const GenomePEDData& ped_data,
                            std::ostream& outfile);

 static bool writeSynResults(const ContigId_t& contig_id,
                             const ResultsMap& genome_results_map,
                             std::ostream& outfile);

};



} // namespace





#endif //KGL_ANALYSIS_MUTATION_INBREED_H
