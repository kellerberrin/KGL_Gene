//
// Created by kellerberrin on 7/12/20.
//

#include "kgl_analysis_mutation_inbreed_execute.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_analysis_mutation_synthetic.h"
#include "kgl_database/kgl_filter.h"


namespace kgl = kellerberrin::genome;



// Perform the genetic analysis per iteration.
bool kgl::ExecuteInbreedingAnalysis::executeAnalysis(std::shared_ptr<const GenomeReference> genome_GRCh38,
                                                     std::shared_ptr<const DiploidPopulation> diploid_population,
                                                     std::shared_ptr<const UnphasedPopulation> unphased_population,
                                                     std::shared_ptr<const GenomePEDData> ped_data) {

  genome_GRCh38_ = std::move(genome_GRCh38);
  diploid_population_ = std::move(diploid_population);
  unphased_population_ = std::move(unphased_population);
  ped_data_ = std::move(ped_data);

  // If necessary, create an unphased population from the diploid population.
  createUnphased();

  if (analyze_diploid_) {

    // Check that we have what we need.
    if (not (diploid_population_ and unphased_population_ and ped_data_)) {

      ExecEnv::log().error("ExecuteInbreedingAnalysis::processDiploid; Insufficient data, cannot process diploid inbreeding");
      return false;

    }

    return processDiploid();

  } else {

    // Check that we have what we need.
    if (not unphased_population_) {

      ExecEnv::log().error("MutationAnalysis::iterationAnalysis; Insufficient data, cannot process synthetic diploid inbreeding");
      return false;

    }

    return processSynthetic();

  }

}

void kgl::ExecuteInbreedingAnalysis::createUnphased() {

  if (diploid_population_ and not unphased_population_) {

    ExecEnv::log().info("ExecuteInbreedingAnalysis::processDiploid; Creating unique unphased population using 1000 Genomes.");
    std::shared_ptr<GenomeVariant> unphased_genome_ptr = diploid_population_->uniqueUnphasedGenome<GenomeVariant>();
    std::shared_ptr<UnphasedPopulation> unphased_unique_ptr = std::make_shared<UnphasedPopulation>(diploid_population_->populationId());
    unphased_unique_ptr->addGenome(unphased_genome_ptr);
    ExecEnv::log().info("ExecuteInbreedingAnalysis::processDiploid; Created unique unphased population, variant count: {}.", unphased_unique_ptr->variantCount());
    unphased_population_ = unphased_unique_ptr;

  }

}


// Perform the genetic analysis per iteration.
bool kgl::ExecuteInbreedingAnalysis::processDiploid() {


  // Setup the analysis parameters.
  size_t locii_count = 1000000000;
  static const ContigOffset_t sampling_distance = 10000;
  static const ContigOffset_t lower_window = 0;
  static const ContigOffset_t upper_window = 1000000000;
  static const double upper_allele_frequency = 1.0;
  static const double lower_allele_frequency = 0.20;

  // Filter out any variants that did not pass VCF filters (otherwise we get duplicate variants).
  unphased_population_ = unphased_population_->filterVariants(AndFilter(SNPFilter(), PassFilter()));

  InbreedingParameters parameters;
  parameters.lociiArguments().minAlleleFrequency(lower_allele_frequency);
  parameters.lociiArguments().maxAlleleFrequency(upper_allele_frequency);
  parameters.inbreedingAlgorthim(InbreedingCalculation::LOGLIKELIHOOD_F);
  parameters.lociiArguments().frequencySource(FrequencyDatabaseRead::alleleFrequencySource(unphased_population_));
  parameters.lociiArguments().lowerOffset(lower_window);
  parameters.lociiArguments().upperOffset(upper_window);
  parameters.lociiArguments().lociiCount(locii_count);
  parameters.lociiArguments().lociiSpacing(sampling_distance);

  InbreedingOutputResults results(diploid_population_->populationId());

  InbreedingAnalysis::populationInbreeding( unphased_population_,
                                            *diploid_population_,
                                            *ped_data_,
                                            parameters,
                                            results);

  inbreeding_results_.emplace(results.identifier(), results);

  return true;

}

// Perform an inbreeding analysis of a synthetic population.
bool kgl::ExecuteInbreedingAnalysis::processSynthetic() {

  // Setup the analysis parameters.
  size_t locii_count = 3000;
  static const ContigOffset_t sampling_distance = 0;
  static const ContigOffset_t lower_window = 115000000;
  static const ContigOffset_t upper_window = 175000000;
  static const double upper_allele_frequency = 1.0;
  static const double lower_allele_frequency = 0.2;

  InbreedingParameters parameters;
  parameters.lociiArguments().minAlleleFrequency(lower_allele_frequency);
  parameters.lociiArguments().maxAlleleFrequency(upper_allele_frequency);
  parameters.inbreedingAlgorthim(InbreedingCalculation::LOGLIKELIHOOD_F);
  parameters.lociiArguments().frequencySource(FrequencyDatabaseRead::alleleFrequencySource(unphased_population_));
  parameters.lociiArguments().lowerOffset(lower_window);
  parameters.lociiArguments().upperOffset(upper_window);
  parameters.lociiArguments().lociiCount(locii_count);
  parameters.lociiArguments().lociiSpacing(sampling_distance);

  InbreedingOutputResults results(unphased_population_->populationId());

  SyntheticAnalysis::syntheticInbreeding(unphased_population_, parameters, results);

  inbreeding_results_.emplace(results.identifier(), results);

  return true;

}


bool kgl::ExecuteInbreedingAnalysis::writeResults(const std::string& output_file) {


  for (auto const& [identifier, result] :  inbreeding_results_) {

    ExecEnv::log().info("ExecuteInbreedingAnalysis::processDiploid: Writing results for identifier: {}", identifier);

    if (analyze_diploid_) {

      InbreedingOutput::writeColumnResults(result, *ped_data_, output_file);

    } else {

      InbreedingOutput::writeSynResults(result, output_file);

    }

  }

  return true;

}

