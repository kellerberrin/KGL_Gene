//
// Created by kellerberrin on 7/12/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_EXECUTE_H
#define KGL_ANALYSIS_MUTATION_INBREED_EXECUTE_H


#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"
#include "kgl_analysis_mutation_inbreed_calc.h"
#include "kgl_analysis_mutation_inbreed_output.h"


namespace kellerberrin::genome {   //  organization::project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Sets up the Inbreeding Parameters and then executes the analysis.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using InbreedingResultMap = std::map<std::string, InbreedingOutputResults>;

class ExecuteInbreedingAnalysis {

public:

  // The population variant data.
  ExecuteInbreedingAnalysis() = default;

  ~ExecuteInbreedingAnalysis() = default;

  bool executeAnalysis(std::shared_ptr<const GenomeReference> genome_GRCh38,
                       std::shared_ptr<const DiploidPopulation> diploid_population,
                       std::shared_ptr<const UnphasedPopulation> unphased_population,
                       std::shared_ptr<const GenomePEDData> ped_data);

  bool writeResults(const std::string &output_file);

private:

  // The population variant data.
  std::shared_ptr<const GenomeReference> genome_GRCh38_;
  std::shared_ptr<const DiploidPopulation> diploid_population_;
  std::shared_ptr<const UnphasedPopulation> unphased_population_;
  std::shared_ptr<const GenomePEDData> ped_data_;
  // Analysis results generated to be output to file on finalize.
  InbreedingResultMap inbreeding_results_;

  const static bool analyze_diploid_ = true;

  bool processDiploid();

  bool processSynthetic();

  void createUnphased();

  // Required because different allele frequency sources have
  // different allele frequency codes for super populations.
  FrequencyDatabaseSource alleleFrequencySource();


};


} // namespace




#endif //KGL_ANALYSIS_MUTATION_INBREED_EXECUTE_H
