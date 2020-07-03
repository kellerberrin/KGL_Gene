//
// Created by kellerberrin on 3/7/20.
//

#include "kgl_analysis_mutation.h"
#include "kgl_variant_db_phased.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::MutationAnalysis::initializeAnalysis(const std::string& work_directory,
                                               const RuntimeParameterMap& named_parameters,
                                               std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Default Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::fileReadAnalysis(std::shared_ptr<const PopulationBase> population_base) {

  // Superclass the population
  std::shared_ptr<const DiploidPopulation> population = std::dynamic_pointer_cast<const DiploidPopulation>(population_base);

  if (not population) {

    ExecEnv::log().error("Analysis: {},  expected a Diploid Population", ident());
    return false;

  }

  ExecEnv::log().info("Analysis: {}, initialized", ident());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::MutationAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}

