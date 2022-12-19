//
// Created by kellerberrin on 26/7/21.
//

#include "kgl_analysis_json.h"
#include "kgl_json_parser.h"

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

  // Truncate the file to zero size.
  std::ofstream citation_file(citation_file_name_);

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::JsonAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("Analysis Id: {}, read file: {}", ident(), data_ptr->fileId());

  auto file_characteristic = data_ptr->dataCharacteristic();

  if (file_characteristic.data_implementation == DataImplEnum::PMIDCitationMap) {

    std::shared_ptr<const DBCitation> citation_ptr = std::dynamic_pointer_cast<const DBCitation>(data_ptr);
    if (not citation_ptr) {

      ExecEnv::log().critical("JsonAnalysis::fileReadAnalysis, provided file type is not a DBCitation, unrecoverable.");

    }

    if (not writeAppendCitations(citation_ptr->citationMap())) {

      ExecEnv::log().critical("JsonAnalysis::fileReadAnalysis; problem processing file: {}", data_ptr->fileId());

    }

  } else if (file_characteristic.data_implementation == DataImplEnum::FileName) {

    ExecEnv::log().info("JsonAnalysis::fileReadAnalysis; Json file name specified : {}", data_ptr->fileId());
    json_file_names_.push_back(data_ptr->fileId());

  } else {

    ExecEnv::log().critical("JsonAnalysis::fileReadAnalysis; Unknown Data Type Passed to Analysis Id: {}", ident());
    return false;

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

  ExecEnv::log().info("Finalize Analysis called for Analysis Id: {}", ident());

  WorkflowThreads thread_pool(WorkflowThreads::defaultThreads(json_file_names_.size()));
  // Simplify the future vector type definition.
  using CitationResult = std::shared_ptr<const DBCitationMap>;
  using FutureResult = std::future<CitationResult>;
  std::vector<FutureResult> future_vector;

  for (auto const& json_file : json_file_names_) {

    FutureResult future = thread_pool.enqueueTask(&JsonAnalysis::parseJsonFile, json_file);

    future_vector.push_back(std::move(future));

  }

  for (auto &future : future_vector) {

    auto citation_map_ptr = future.get();
    if (not writeAppendCitations(*citation_map_ptr)) {

      ExecEnv::log().error("JsonAnalysis::finalizeAnalysis; problem writing to Citation File");

    }

  }

  return true;

}


std::shared_ptr<const kgl::DBCitationMap> kgl::JsonAnalysis::parseJsonFile(std::string json_file) {

  std::shared_ptr<DBCitationMap> citation_map_ptr(std::make_shared<DBCitationMap>());

  JSONInfoParser reader;
  if (not reader.parseFile(json_file, *citation_map_ptr)) {

    ExecEnv::log().error("JsonAnalysis::parseJsonFile, error parsing Json file: {}", json_file);

  }

  return citation_map_ptr;

}


bool kgl::JsonAnalysis::writeAppendCitations(const DBCitationMap& citation_map) const {

  // Append to existing file.
  std::ofstream citation_file(citation_file_name_ , std::fstream::out | std::fstream::app);

  if (not citation_file.good()) {

    ExecEnv::log().error("kgl::JsonAnalysis::writeAppendCitations, error opening citation file: {}", citation_file_name_);
    return false;

  }

  // Append the citation map to file.
  size_t citation_count{0};
  for (auto const& [rsid, citations] : citation_map) {

    for (auto const& citation : citations) {

      ++citation_count;
      citation_file << rsid << OUTPUT_DELIMITER_ << citation << '\n';

    }

  }

  ExecEnv::log().info("Successfully Appended alleles: {}, citations: {} to file: {}", citation_map.size(), citation_count, citation_file_name_);

  return true;

}