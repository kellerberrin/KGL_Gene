//
// Created by kellerberrin on 3/7/20.
//

#ifndef KGL_ANALYSIS_MUTATION_H
#define KGL_ANALYSIS_MUTATION_H


#include "kgl_analysis_virtual.h"
#include "kgl_ped_parser.h"
#include "kgl_analysis_mutation_inbreed_execute.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class MutationAnalysis : public VirtualAnalysis {

public:

  MutationAnalysis() = default;

  ~MutationAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return ANALYSIS_IDENT_; }

  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<MutationAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis(const std::string &work_directory,
                                        const RuntimeParameterMap &named_parameters,
                                        std::shared_ptr<const GenomeCollection> reference_genomes) override;

  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  constexpr static const char* REFERENCE_GENOME_ = "GRCh38";
  constexpr static const char* OUTPUT_FILE_ = "OUTPUTFILE";
  constexpr static const char* ANALYSIS_IDENT_ = "MUTATION";

  [[nodiscard]] bool getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters);
  std::string output_file_name_;

  // The population variant data.
  std::shared_ptr<const GenomeReference> genome_GRCh38_;
  std::shared_ptr<const PopulationDB> diploid_population_;
  std::shared_ptr<const PopulationDB> unphased_population_;
  std::shared_ptr<const GenomePEDData> ped_data_;

  ExecuteInbreedingAnalysis inbreed_analysis_;

};



} // namespace

#endif //KGL_KGL_ANALYSIS_MUTATION_H
