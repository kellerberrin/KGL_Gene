//
// Created by kellerberrin on 3/1/21.
//

#ifndef KGL_ANALYSIS_PFEMP_H
#define KGL_ANALYSIS_PFEMP_H


#include "kgl_analysis_virtual.h"



namespace kellerberrin::genome {   //  organization::project level namespace


class PfEMPAnalysis : public VirtualAnalysis {

public:

  PfEMPAnalysis() = default;
  ~PfEMPAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "PfEMP"; }
  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<PfEMPAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const ActiveParameterList& named_parameters,
                                         std::shared_ptr<const GenomeCollection> reference_genomes) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  std::shared_ptr<const GenomeCollection> reference_genomes_;
  std::string work_directory_;

  void performPFEMP1UPGMA();


};


} // namespace








#endif //KGL_KGL_ANALYSIS_PFEMP_H
