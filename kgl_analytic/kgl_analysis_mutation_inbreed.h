//
// Created by kellerberrin on 23/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_H
#define KGL_ANALYSIS_MUTATION_INBREED_H


#include "kgl_analysis_virtual.h"
#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class InbreedingAnalysis  {

public:

  InbreedingAnalysis() = delete;
  ~InbreedingAnalysis() = delete;

  // Construct a synthetic population and analyze it.
  [[nodiscard]] static bool syntheticInbreeding( const GenomeReference& genome_GRCh38,
                                                 const UnphasedPopulation& unphased_population,
                                                 const std::string& output_file_name);

  // Process Inbreeding coefficient and relatedness using pre-generated allele locus lists.
  [[nodiscard]] static bool populationInbreeding(const GenomeReference& genome_GRCh38,
                                                 const UnphasedPopulation& unphased_population,
                                                 const DiploidPopulation& diploid_population,
                                                 const GenomePEDData& ped_data,
                                                 const std::string& output_file_name);

private:

  // Synthetic inbreeding values.
  constexpr static const double MIN_INBREEDING_COEFICIENT = -0.5; // Set the minimum inbreeding coefficient
  constexpr static const double MAX_INBREEDING_COEFICIENT = 0.5; // Set the maximum inbreeding coefficient
  constexpr static const double STEP_INBREEDING_COEFICIENT = 0.01; // Set the inbreeding step


  // Called by the threadpool for each genome/sample.
  struct LocusResults {

    GenomeId_t genome;
    size_t hetero_count{0};
    size_t homo_count{0};
    size_t total_allele_count{0};
    double inbred_allele_sum{0.0};

  };
  using ResultsMap = std::map<GenomeId_t, LocusResults>;

  // The inbreeding Algorithms
  [[nodiscard]] static LocusResults processRitlandLocus( const ContigId_t& genome_id,
                                                         const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                         const std::string& super_population_field,
                                                         const std::shared_ptr<const ContigVariant>& locus_list);

  [[nodiscard]] static LocusResults multiLocus1( const ContigId_t& contig_id,
                                                 const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                 const std::string& super_population_field,
                                                 const std::shared_ptr<const ContigVariant>& locus_list);

  constexpr static const char DELIMITER_ = ',';
  constexpr static const char* FILE_EXT_ = ".csv";

  // Write the analysis results to a CSV file.
  static bool writeResults( const ContigId_t& contig_id,
                            const ResultsMap& locus_results,
                            const GenomePEDData& ped_data,
                            const std::string& output_file_name);

 static bool syntheticResults( const ContigId_t& contig_id,
                               const ResultsMap& genome_results_map,
                               const std::string& output_file_name);

};



} // namespace





#endif //KGL_ANALYSIS_MUTATION_INBREED_H
