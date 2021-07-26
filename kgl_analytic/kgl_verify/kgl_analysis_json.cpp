//
// Created by kellerberrin on 26/7/21.
//

#include "kgl_analysis_json.h"

namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::JsonAnalysis::initializeAnalysis( const std::string& work_directory,
                                              const ActiveParameterList& named_parameters,
                                              const std::shared_ptr<const AnalysisResources>&) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident);

  }

  citation_file_name_ = std::string(CITATION_OUTPUT_FILE_) + std::string(OUTPUT_FILE_EXT_);
  citation_file_name_ = Utility::filePath(citation_file_name_ , work_directory);

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::JsonAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {}", ident());

  auto file_characteristic = data_ptr->dataCharacteristic();

  if (file_characteristic.data_implementation == DataImplEnum::PMIDCitationMap) {

    std::shared_ptr<const DBCitation> citation_ptr = std::dynamic_pointer_cast<const DBCitation>(data_ptr);
    if (not citation_ptr) {

      ExecEnv::log().critical("JsonAnalysis::fileReadAnalysis, provided file type is not a DBCitation, unrecoverable.");

    }

    // Add the citations to the analysis map
    for (auto const& [rsid, citations] : citation_ptr->citationMap()) {

      auto [iter, result] = citation_map_.try_emplace(rsid, citations);

      if (not result) {

        ExecEnv::log().warn("JsonAnalysis::fileReadAnalysis; unable to add (duplicate) rsid: {} to citation map", rsid);

      }

    }

  } else {


    ExecEnv::log().warn("JsonAnalysis::fileReadAnalysis; unexpected file: {}", data_ptr->fileId());

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::JsonAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All Json files have been presented, write citation file as a 2 column flat file.
bool kgl::JsonAnalysis::finalizeAnalysis() {

  std::ofstream citation_file(citation_file_name_);

  if (not citation_file.good()) {

    ExecEnv::log().error("JsonAnalysis::finalizeAnalysis, error openning citation file: {}", citation_file_name_);
    return false;

  }

  // Write the citation map to file.
  size_t citation_count{0};
  for (auto const& [rsid, citations] : citation_map_) {

    for (auto const& citation : citations) {

      ++citation_count;
      citation_file << rsid << OUTPUT_DELIMITER_ << citation << '\n';

    }

  }

  ExecEnv::log().info("Successfully wrote alleles: {}, citations: {} to file: {}", citation_map_.size(), citation_count, citation_file_name_);

  return true;

}

