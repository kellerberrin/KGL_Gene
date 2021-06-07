//
// Created by kellerberrin on 3/1/21.
//

#ifndef KGL_ANALYSIS_PFEMP_H
#define KGL_ANALYSIS_PFEMP_H


#include "kgl_analysis_virtual.h"
#include "kgl_genome_collection.h"



namespace kellerberrin::genome {   //  organization::project level namespace


class PfEMPAnalysis : public VirtualAnalysis {

public:

  PfEMPAnalysis() { reference_genomes_ = std::make_shared<GenomeCollection>(); }
  ~PfEMPAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "PfEMP"; }
  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<PfEMPAnalysis>(); }

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

  std::shared_ptr<GenomeCollection> reference_genomes_;

  void performPFEMP1UPGMA();
  bool getParameters(const ActiveParameterList& named_parameters, const std::string& work_directory);

  // Available Parameters
  constexpr static const char* NEWICK_FILE_ = "NewickFile";
  constexpr static const char* INTRON_FILE_ = "IntronFile";

  std::string newick_file_name_;
  std::string intron_file_name_;

};


} // namespace








#endif //KGL_KGL_ANALYSIS_PFEMP_H
