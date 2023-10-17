//
// Created by kellerberrin on 4/5/20.
//

#ifndef KGL_NULL_ANALYSIS_H
#define KGL_NULL_ANALYSIS_H


#include "kgl_package_analysis_virtual.h"
#include "kgl_pubmed_resource.h"


namespace kellerberrin::genome::analysis {   //  organization::project level namespace


class NullAnalysis : public VirtualAnalysis {

public:

  NullAnalysis() = default;
  ~NullAnalysis() override = default;

  // The ident must match the ident used in the package XML.
  constexpr static std::string IDENT {"NULL"};
  // Need a polymorphic version to interrogate VirtualAnalysis pointers.
  [[nodiscard]] std::string ident() const override { return IDENT; }
  // Simple creation factory function.
  [[nodiscard]] static std::unique_ptr<VirtualAnalysis> factory() { return std::make_unique<NullAnalysis>(); }


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


};


} // namespace

#endif //KGL_NULL_ANALYSIS_H
