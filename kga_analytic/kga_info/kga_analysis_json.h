//
// Created by kellerberrin on 26/7/21.
//

#ifndef KGL_ANALYSIS_JSON_H
#define KGL_ANALYSIS_JSON_H


#include "kgl_package_analysis_virtual.h"
#include "kgl_json_parser.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object parses Allele Json Files and typically writes the parsed information to a resource file.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome::analysis {   //  organization::project level namespace


class JsonAnalysis : public VirtualAnalysis {

public:

  JsonAnalysis() = default;
  ~JsonAnalysis() override = default;

  // The ident must match the ident used in the package XML.
  constexpr static std::string IDENT {"PARSEJSON"};
  // Need a polymorphic version to interrogate VirtualAnalysis pointers.
  [[nodiscard]] std::string ident() const override { return IDENT; }
  // Simple creation factory function.
  [[nodiscard]] static std::unique_ptr<VirtualAnalysis> factory() { return std::make_unique<JsonAnalysis>(); }

  // Setup the analytics to process data.
  [[nodiscard]] bool initializeAnalysis(const std::string &work_directory,
                                        const ActiveParameterList &named_parameters,
                                        const std::shared_ptr<const AnalysisResources> &resource_ptr) override;

  // Perform the genetic analysis per file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  std::string citation_file_name_;
  std::vector<std::string> json_file_names_;

  constexpr static const char* CITATION_OUTPUT_FILE_ = "CitationMT";
  constexpr static const char OUTPUT_DELIMITER_ = ',';
  constexpr static const char* OUTPUT_FILE_EXT_ = ".csv";

  [[nodiscard]] bool writeAppendCitations(const DBCitationMap& citation_map) const;
  [[nodiscard]] static std::shared_ptr<const DBCitationMap> parseJsonFile(std::string json_file);



  };


} // namespace


#endif //KGL_ANALYSIS_JSON_H
