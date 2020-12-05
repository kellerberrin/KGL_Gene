//
// Created by kellerberrin on 3/7/20.
//

#include "kel_thread_pool.h"
#include "kel_distribution.h"
#include "kgl_analysis_mutation.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_analysis_mutation_synthetic.h"
#include "kel_optimize.h"

#include <fstream>

namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::MutationAnalysis::initializeAnalysis(const std::string& work_directory,
                                               const RuntimeParameterMap& named_parameters,
                                               std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  std::optional<std::shared_ptr<const GenomeReference>> ref_genome_opt = reference_genomes->getOptionalGenome(REFERENCE_GENOME_);

  if (ref_genome_opt) {

    genome_GRCh38_ = ref_genome_opt.value();

  } else {

    ExecEnv::log().error("MutationAnalysis::initializeAnalysis, Could not find Genome: {} Analysis: {} disabled.", REFERENCE_GENOME_, ident());
    return false;

  }

  if (not getParameters(work_directory, named_parameters)) {

    return false;

  }

  return true;

}

// This function superclasses the data objects and stores them for further use.
bool kgl::MutationAnalysis::fileReadAnalysis(std::shared_ptr<const DataObjectBase> data_object_ptr) {

  ExecEnv::log().info("Analysis: {}, begin processing data file", ident(), data_object_ptr->Id());


  if (data_object_ptr->dataType() == DataTypeEnum::DiploidPopulation) {

    std::shared_ptr<const DiploidPopulation> diploid_population = std::dynamic_pointer_cast<const DiploidPopulation>(data_object_ptr);

    if (diploid_population) {

      ExecEnv::log().info("Analysis: {}, Generate inbreeding statistics for file: {}", ident(), data_object_ptr->Id());

      diploid_population_ = diploid_population;

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a Diploid Population", ident(), data_object_ptr->Id());
      return false;

    }

  }

  if (data_object_ptr->dataType() == DataTypeEnum::UnphasedPopulation) {

    std::shared_ptr<const UnphasedPopulation> unphased_population = std::dynamic_pointer_cast<const UnphasedPopulation>(data_object_ptr);

    if (unphased_population) {

      unphased_population_ = unphased_population;

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not an Unphased Population", ident(), data_object_ptr->Id());
      return false;

    }

  }

  if (data_object_ptr->dataType() == DataTypeEnum::PedAncestor) {

    std::shared_ptr<const GenomePEDData> ped_data = std::dynamic_pointer_cast<const GenomePEDData>(data_object_ptr);

    if (ped_data) {

      ped_data_ = ped_data;
      ExecEnv::log().info("Analysis: {}, ped file: {} contains: {} PED records", ident(), ped_data->Id(), ped_data->getMap().size());

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a PED Ancestor Object", ident(), data_object_ptr->Id());
      return false;

    }

  }

  ExecEnv::log().info("Analysis: {}, completed data file: {}", ident(), data_object_ptr->Id());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Iteration Analysis called for Analysis Id: {}", ident());

  // If necessary, create an unphased population from the diploid population.
  createUnphased();

  if (analyze_diploid_) {

    // Check that we have what we need.
    if (not (diploid_population_ and unphased_population_ and ped_data_)) {

      ExecEnv::log().error("MutationAnalysis::iterationAnalysis; Insufficient data, cannot process diploid inbreeding");
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

void kgl::MutationAnalysis::createUnphased() {

  if (diploid_population_ and not unphased_population_) {

    ExecEnv::log().info("MutationAnalysis::createUnphased; Creating unique unphased population using 1000 Genomes.");
    std::shared_ptr<GenomeVariant> unphased_genome_ptr = diploid_population_->uniqueUnphasedGenome<GenomeVariant>();
    std::shared_ptr<UnphasedPopulation> unphased_unique_ptr = std::make_shared<UnphasedPopulation>(diploid_population_->populationId());
    unphased_unique_ptr->addGenome(unphased_genome_ptr);
    ExecEnv::log().info("MutationAnalysis::createUnphased; Created unique unphased population, variant count: {}.", unphased_unique_ptr->variantCount());
    unphased_population_ = unphased_unique_ptr;

  }

}

kgl::FrequencyDatabaseSource kgl::MutationAnalysis::alleleFrequencySource() {

  static const std::string gnomad_3_1_fragment = "r3.1";
  static const std::string gnomad_2_1_fragment = "r2.1";
  static const std::string genome_1000_fragment = "1000";

  size_t find_pos = unphased_population_->populationId().find(genome_1000_fragment);
  if (find_pos != std::string::npos) {

    return FrequencyDatabaseSource::GENOMES_1000;

  }

  find_pos = unphased_population_->populationId().find(gnomad_2_1_fragment);
  if (find_pos != std::string::npos) {

    return FrequencyDatabaseSource::GNOMAD2_1;

  }

  find_pos = unphased_population_->populationId().find(gnomad_3_1_fragment);
  if (find_pos != std::string::npos) {

    return FrequencyDatabaseSource::GNOMAD3_1;

  }

  // Id signature not found, complain and return GNOMAD2_1.
  ExecEnv::log().error("MutationAnalysis::alleleFrequencySource; Population allele frequency signature not for unphased population: {}",
                       unphased_population_->populationId());

  return FrequencyDatabaseSource::GNOMAD2_1;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::processDiploid() {


  // Setup the analysis parameters.
  size_t locii_count = 3000;
  static const ContigOffset_t sampling_distance = 100;
  static const ContigOffset_t lower_window = 0;
  static const ContigOffset_t upper_window = 1000000000;
  static const double upper_allele_frequency = 1.0;
  static const double lower_allele_frequency = 0.2;

  InbreedingParameters parameters;
  parameters.lociiArguments().minAlleleFrequency(lower_allele_frequency);
  parameters.lociiArguments().maxAlleleFrequency(upper_allele_frequency);
  parameters.inbreedingAlgorthim(InbreedingCalculation::LOGLIKELIHOOD);
  parameters.lociiArguments().frequencySource(alleleFrequencySource());
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
bool kgl::MutationAnalysis::processSynthetic() {

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
  parameters.inbreedingAlgorthim(InbreedingCalculation::LOGLIKELIHOOD);
  parameters.lociiArguments().frequencySource(alleleFrequencySource());
  parameters.lociiArguments().lowerOffset(lower_window);
  parameters.lociiArguments().upperOffset(upper_window);
  parameters.lociiArguments().lociiCount(locii_count);
  parameters.lociiArguments().lociiSpacing(sampling_distance);

  InbreedingOutputResults results(unphased_population_->populationId());

  SyntheticAnalysis::syntheticInbreeding(unphased_population_, parameters, results);

  inbreeding_results_.emplace(results.identifier(), results);

  return true;

}


// All VCF data has been presented, finalize analysis and write results.
bool kgl::MutationAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Finalize called for Analysis Id: {}", ident());

  for (auto const& [identifier, result] :  inbreeding_results_) {

    ExecEnv::log().info("Analysis Id: {}; Writing results for identifier: {}", ident(), identifier);

    if (analyze_diploid_) {

      InbreedingOutput::writeColumnResults(result, *ped_data_, output_file_name_);

    } else {

      InbreedingOutput::writeSynResults(result, output_file_name_);

    }

  }

  return true;

}


bool kgl::MutationAnalysis::getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters) {

  // Get the output filename
  auto result = named_parameters.find(OUTPUT_FILE_);
  if (result == named_parameters.end()) {
    ExecEnv::log().error("Analytic: {}; Expected Parameter: {} to be defined. {} is deactivated. Available named Parameters:", ident(), OUTPUT_FILE_, ident());
    for (auto const& [parameter_ident, parameter_value] : named_parameters) {

      ExecEnv::log().info("Analysis: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

    }
    return false;
  }
  output_file_name_ = Utility::filePath(result->second, work_directory);

  ExecEnv::log().info("Analysis: {}, initialized with output file: {}", ident(), output_file_name_);

  return true;

}



