//
// Created by kellerberrin on 5/1/21.
//

#ifndef KGL_ANALYSIS_GENE_H
#define KGL_ANALYSIS_GENE_H


#include "kgl_analysis_virtual.h"
#include "kgl_Pf3k_COI.h"



namespace kellerberrin::genome {   //  organization::project level namespace


class MutationAnalysis : public VirtualAnalysis {

public:

  MutationAnalysis() = default;
  ~MutationAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "MUTATION"; }
  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<MutationAnalysis>(); }

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


  std::shared_ptr<const GenomeCollection> genome_collection_ptr_;
  std::shared_ptr<const PopulationDB> population_ptr_;
  std::shared_ptr<const Pf3kCOIDB> pf3k_coi_ptr_;

  void performRegion();

  constexpr static const char* GENOME_IDENT_ = "AnalysisGenome";

  constexpr static const char* OUTPUT_FILE_ = "OutputFile";
  constexpr static const char OUTPUT_DELIMITER_ = ',';
  constexpr static const char* OUTPUT_FILE_EXT_ = ".csv";

  std::string work_directory_;
  std::string output_file_name_;

  bool getParameters(const ActiveParameterList& named_parameters,
                     const std::shared_ptr<const GenomeCollection>& genome_collection_ptr,
                     const std::string& work_directory);
  bool writeOutput();
  void writeHeader(std::ostream& out_file);

};


} // namespace


#endif //KGL_ANALYSIS_GENE_H
