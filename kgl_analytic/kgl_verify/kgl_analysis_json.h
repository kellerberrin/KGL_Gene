//
// Created by kellerberrin on 26/7/21.
//

#ifndef KGL_ANALYSIS_JSON_H
#define KGL_ANALYSIS_JSON_H


#include "kgl_analysis_virtual.h"
#include "kgl_Json_parser.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object parses Allele Json Files and typically writes the information to file.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization::project level namespace


class JsonAnalysis : public VirtualAnalysis {

public:

  JsonAnalysis() = default;
  ~JsonAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "PARSEJSON"; }

  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<JsonAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis(const std::string &work_directory,
                                        const ActiveParameterList &named_parameters,
                                        const std::shared_ptr<const AnalysisResources> &resource_ptr) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  DBCitationMap citation_map_;
  std::string citation_file_name_;

  constexpr static const char* CITATION_OUTPUT_FILE_ = "CitationOut";
  constexpr static const char OUTPUT_DELIMITER_ = ',';
  constexpr static const char* OUTPUT_FILE_EXT_ = ".csv";

};


} // namespace


#endif //KGL_ANALYSIS_JSON_H
