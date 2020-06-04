//
// Created by kellerberrin on 4/5/20.
//

#include "kgl_analysis_null.h"

namespace kgl = kellerberrin::genome;

// Setup the analytics to process VCF data.
bool kgl::NullAnalysis::initializeAnalysis( const std::string& work_directory,
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
bool kgl::NullAnalysis::fileReadAnalysis(std::shared_ptr<const UnphasedPopulation> population) {

  ExecEnv::log().info("Default VCF File Read for Analysis Id: {} called with Variant Population: {}, Variant Count: {}", ident(), population->populationId(), population->variantCount());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::NullAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::NullAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}
