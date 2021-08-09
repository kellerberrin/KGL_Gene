//
// Created by kellerberrin on 4/5/20.
//

#ifndef KGL_NULL_ANALYSIS_H
#define KGL_NULL_ANALYSIS_H


#include "kgl_analysis_virtual.h"
#include "kel_rest_api.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class NullAnalysis : public VirtualAnalysis {

public:

  NullAnalysis() = default;
  ~NullAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "NULL"; }
  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<NullAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const ActiveParameterList& named_parameters,
                                         const std::shared_ptr<const AnalysisResources>& resource_ptr) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  std::string work_directory_;
  RestAPI test_rest_api_;

  void investigateVepFields(const std::shared_ptr<const PopulationDB>& population);
  void genomeIdIndex(const std::shared_ptr<const PopulationDB>& population);
  void ensemblIdIndex(const std::shared_ptr<const PopulationDB>& population);

};


} // namespace

#endif //KGL_NULL_ANALYSIS_H
