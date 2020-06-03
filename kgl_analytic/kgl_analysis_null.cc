//
// Created by kellerberrin on 4/5/20.
//

#include "kgl_analysis_null.h"

namespace kgl = kellerberrin::genome;

// Setup the analytics to process VCF data.
bool kgl::NullAnalysis::initializeAnalysis( const std::string& work_directory,
                                            const RuntimeParameterMap& named_parameters,
                                            std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Analysis: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters) {

    ExecEnv::log().info("Initialize Analysis: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Initialize Analysis {} called with Reference Genome: {}", ident(), genome.first);

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::NullAnalysis::fileReadAnalysis(std::shared_ptr<const GenomeCollection>  reference_genomes, std::shared_ptr<const UnphasedPopulation> population) {


  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("File Read Analysis: {} called with Reference Genome: {}", ident(), genome.first);

  }

  ExecEnv::log().info("File Read Analysis: {} called with Variant Population: {}, Variant Count: {}", ident(), population->populationId(), population->variantCount());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::NullAnalysis::iterationAnalysis(std::shared_ptr<const GenomeCollection>  reference_genomes) {


  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Iteration Analysis: {} called with Reference Genome: {}", ident(), genome.first);

  }

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::NullAnalysis::finalizeAnalysis(std::shared_ptr<const GenomeCollection> reference_genomes) {

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Finalize Analysis: {} called with Reference Genome: {}", ident(), genome.first);

  }

  return true;

}
