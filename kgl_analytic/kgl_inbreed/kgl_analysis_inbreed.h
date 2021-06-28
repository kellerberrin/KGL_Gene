//
// Created by kellerberrin on 3/7/20.
//

#ifndef KGL_ANALYSIS_MUTATION_H
#define KGL_ANALYSIS_MUTATION_H


#include "kgl_analysis_virtual.h"
#include "kgl_genealogy_parser.h"
#include "kgl_analysis_inbreed_args.h"
#include "kgl_analysis_inbreed_output.h"
#include "kgl_analysis_inbreed_execute.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class InbreedAnalysis : public VirtualAnalysis {

public:

  InbreedAnalysis() = default;

  ~InbreedAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return ANALYSIS_IDENT_; }

  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<InbreedAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis(const std::string &work_directory,
                                        const ActiveParameterList& named_parameters,
                                        const std::shared_ptr<const AnalysisResources>& resource_ptr) override;

  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  constexpr static const char* ANALYSIS_IDENT_ = "INBREED";

  // The analysis will be executed for each parameter block.
  std::vector<InbreedParamOutput> parameter_output_vector_;
  std::string work_directory_;
  // The population variant data.
  std::shared_ptr<const PopulationDB> diploid_population_;
  std::shared_ptr<const PopulationDB> unphased_population_;
  std::shared_ptr<const GenomeGenealogyData> genealogy_data_;

  // Write to output files.
  bool writeResults();

};



} // namespace

#endif //KGL_KGL_ANALYSIS_MUTATION_H
