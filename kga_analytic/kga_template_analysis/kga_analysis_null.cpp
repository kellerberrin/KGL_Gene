//
// Created by kellerberrin on 4/5/20.
//

#include "kga_analysis_null.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::NullAnalysis::initializeAnalysis(const std::string& work_directory,
                                           const ActiveParameterList& named_parameters,
                                           const std::shared_ptr<const AnalysisResources>& analysis_resource_ptr) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  for (auto const& resource_ptr : analysis_resource_ptr->getResources(ResourceProperties::GENOME_RESOURCE_ID_)) {

    ExecEnv::log().info("Initialize for Analysis Id: {} called with Resource: {}", resource_ptr->resourceIdent());

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::NullAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called with file: {}", ident(), data_ptr->fileId());

  auto file_characteristic = data_ptr->dataCharacteristic();


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

