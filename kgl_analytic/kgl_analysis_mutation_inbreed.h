//
// Created by kellerberrin on 23/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_H
#define KGL_ANALYSIS_MUTATION_INBREED_H


#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"
#include "kgl_analysis_mutation_inbreed_calc.h"
#include "kgl_analysis_mutation_inbreed_locus.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class InbreedingAnalysis  {

public:

  InbreedingAnalysis() = delete;
  ~InbreedingAnalysis() = delete;

  // Analyze population and synthetic data.
  [[nodiscard]] static bool Inbreeding( InbreedingAlgorithm algorithm,
                                        std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                        std::shared_ptr<const DiploidPopulation> diploid_population,
                                        std::shared_ptr<const GenomePEDData> ped_data,
                                        const std::string& output_file_name);

  // Construct a synthetic population and analyze it.
  [[nodiscard]] static bool syntheticInbreeding( InbreedingAlgorithm algorithm,
                                                 std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                 const std::string& output_file_name);

  // Process Inbreeding coefficient and relatedness using pre-generated allele locus lists.
  [[nodiscard]] static bool populationInbreeding(InbreedingAlgorithm algorithm,
                                                 std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                 const DiploidPopulation& diploid_population,
                                                 const GenomePEDData& ped_data,
                                                 const std::string& output_file_name);

private:

  // Synthetic inbreeding values.
  constexpr static const double MIN_INBREEDING_COEFICIENT = -0.5; // Set the minimum inbreeding coefficient
  constexpr static const double MAX_INBREEDING_COEFICIENT = 0.5; // Set the maximum inbreeding coefficient
  constexpr static const double STEP_INBREEDING_COEFICIENT = 0.01; // Set the inbreeding step


  constexpr static const char DELIMITER_ = ',';
  constexpr static const char* FILE_EXT_ = ".csv";

  using ResultsMap = std::map<GenomeId_t, LocusResults>;

  bool static syntheticInbreedingSample( InbreedingAlgorithm algorithm,
                                         std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                         const std::string& output_file_name,
                                         const std::string& sample_name,
                                         double allele_frequency_min,
                                         double allele_frequency_max,
                                         ContigOffset_t  spacing);

  static bool populationInbreedingSample( InbreedingAlgorithm algorithm,
                                          std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                          const DiploidPopulation& diploid_population,
                                          const GenomePEDData& ped_data,
                                          const std::string& output_file_name,
                                          const std::string& sample_name,
                                          double allele_frequency_min,
                                          double allele_frequency_max,
                                          ContigOffset_t  spacing);

  // Use a threadpool to calculate the inbreeding coefficients.
  static bool processResults( InbreedingAlgorithm algorithm,
                              const ContigLocusMap& contig_locus_map,
                              const DiploidPopulation& diploid_population,
                              const GenomePEDData& ped_data,
                              std::ostream& outfile);

  static bool processSynResults( InbreedingAlgorithm algorithm,
                                 const ContigLocusMap& contig_locus_map,
                                 std::ostream& outfile);

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
