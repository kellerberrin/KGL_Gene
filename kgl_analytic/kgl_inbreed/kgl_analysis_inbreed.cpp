//
// Created by kellerberrin on 3/7/20.
//

#include "kel_distribution.h"
#include "kgl_analysis_inbreed.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kel_optimize.h"

#include <fstream>

namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::InbreedAnalysis::initializeAnalysis(const std::string& work_directory,
                                              const ActiveParameterList& named_parameters,
                                              std::shared_ptr<const GenomeCollection>) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters.getMap()) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}", ident(), parameter_ident);

  }

  work_directory_ = work_directory;
  for (auto const& parameter : InbreedArguments::extractParameters(named_parameters)) {

    parameter_output_vector_.emplace_back(InbreedParamOutput(parameter));

  }

  return true;

}

// This function superclasses the data objects and stores them for further use.
bool kgl::InbreedAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) {

  ExecEnv::log().info("Analysis: {}, begin processing data file", ident(), data_object_ptr->fileId());

  auto file_characteristic = data_object_ptr->dataCharacteristic();

  if (file_characteristic.data_structure == DataStructureEnum::DiploidPhased
      or file_characteristic.data_structure == DataStructureEnum::DiploidUnphased) {

    diploid_population_ = std::dynamic_pointer_cast<const PopulationDB>(data_object_ptr);

    if (diploid_population_) {

      ExecEnv::log().info("Analysis: {}, Generate inbreeding statistics for file: {}", ident(), diploid_population_->populationId());

    } else {

      ExecEnv::log().error("InbreedAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a Diploid Population", ident(), data_object_ptr->fileId());
      return false;

    }

  }

  if (file_characteristic.data_structure == DataStructureEnum::UnphasedMonoGenome) {

    unphased_population_  = std::dynamic_pointer_cast<const PopulationDB>(data_object_ptr);

    if (not unphased_population_) {

      ExecEnv::log().error("InbreedAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not an Unphased Population", ident(), data_object_ptr->fileId());
      return false;

    }

    // Only want SNP variants and variants that passed all VCF filters.
    unphased_population_ = unphased_population_->filterVariants(AndFilter(SNPFilter(), PassFilter()));

  }

  if (file_characteristic.data_structure == DataStructureEnum::PedGenome1000) {

    std::shared_ptr<const GenomePEDData> ped_data = std::dynamic_pointer_cast<const GenomePEDData>(data_object_ptr);

    if (ped_data) {

      ped_data_ = ped_data;
      ExecEnv::log().info("Analysis: {}, ped file: {} contains: {} PED records", ident(), ped_data->fileId(), ped_data->getMap().size());

    } else {

      ExecEnv::log().error("InbreedAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a PED Ancestor Object", ident(), data_object_ptr->fileId());
      return false;

    }

  }

  ExecEnv::log().info("Analysis: {}, completed data file: {}", ident(), data_object_ptr->fileId());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::InbreedAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Iteration Analysis called for Analysis Id: {}", ident());

  // If necessary, create an unphased population from the diploid population.
  if (diploid_population_ and not unphased_population_) {

    unphased_population_ = createUnphased();

  }

  for (auto& param_output : parameter_output_vector_) {

    // Set the allele frequency source for this population.
    param_output.getParameters().lociiArguments().frequencySource(unphased_population_->dataSource());
    ExecuteInbreedingAnalysis::executeAnalysis(diploid_population_, unphased_population_, ped_data_, param_output);

  }

  return true;

}



// All VCF data has been presented, finalize analysis and write results.
bool kgl::InbreedAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Finalize called for Analysis Id: {}", ident());

  return writeResults();

}


std::shared_ptr<const kgl::PopulationDB> kgl::InbreedAnalysis::createUnphased() {

  ExecEnv::log().info("InbreedAnalysis::processDiploid; Creating unique unphased population using Diploid Population.");
  std::shared_ptr<GenomeDB> unphased_genome_ptr = diploid_population_->uniqueUnphasedGenome();
  std::shared_ptr<PopulationDB> unphased_unique_ptr = std::make_shared<PopulationDB>(diploid_population_->populationId(),
                                                                                     diploid_population_->dataSource());
  unphased_unique_ptr->addGenome(unphased_genome_ptr);
  ExecEnv::log().info("InbreedAnalysis::processDiploid; Created unique unphased population, variant count: {}.", unphased_unique_ptr->variantCount());

  return unphased_unique_ptr;

}


bool kgl::InbreedAnalysis::writeResults() {


  for (auto& param_output : parameter_output_vector_) {

    if (ped_data_ and not param_output.getParameters().analyzeSynthetic()) {

      InbreedingOutput::writePedResults(param_output, *ped_data_, work_directory_);

    } else {

      InbreedingOutput::writeNoPedResults( param_output, work_directory_);

    }

  }

  return true;

}


