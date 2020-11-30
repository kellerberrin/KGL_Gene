//
// Created by kellerberrin on 3/7/20.
//

#include "kel_thread_pool.h"
#include "kel_distribution.h"
#include "kgl_analysis_mutation.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_analysis_mutation_inbreed.h"
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

  InbreedingParameters parameters;
  parameters.lociiArguments().minAlleleFrequency(0.2);
  parameters.lociiArguments().maxAlleleFrequency(1.0);

  parameters.inbreedingAlgorthim(InbreedingCalculation::LOGLIKELIHOOD);

  InbreedingOutputResults results(diploid_population_->populationId());

  InbreedingAnalysis::InbreedingAll( genome_GRCh38_,
                                     unphased_population_,
                                     diploid_population_,
                                     ped_data_,
                                     parameters,
                                     results);

  InbreedingOutput::writeColumnResults(results , *ped_data_, output_file_name_);

//  inbreeding_results_.emplace(results.identifier(), results);

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::MutationAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Finalize called for Analysis Id: {}", ident());

  for (auto const& [identifier, result] :  inbreeding_results_) {

    ExecEnv::log().info("Analysis Id: {}; Writing results for identifier: {}", ident(), identifier);
    InbreedingOutput::writeColumnResults(result , *ped_data_, output_file_name_);

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



