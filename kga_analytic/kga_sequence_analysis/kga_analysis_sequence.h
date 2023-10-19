//
// Created by kellerberrin on 7/11/20.
//

#ifndef KGL_ANALYSIS_CHECK_H
#define KGL_ANALYSIS_CHECK_H


#include "kgl_package_analysis_virtual.h"
#include "kgl_genome_collection.h"
#include "kga_analysis_sequence_mutation.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object checks the underlying population structures for correctness.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome::analysis {   //  organization::project level namespace


class SequenceAnalysis : public VirtualAnalysis {

public:

  SequenceAnalysis() = default;
  ~SequenceAnalysis() override = default;

  // The ident must match the ident used in the package XML.
  constexpr static std::string IDENT {"VERIFY"};
  // Need a polymorphic version to interrogate VirtualAnalysis pointers.
  [[nodiscard]] std::string ident() const override { return IDENT; }
  // Simple creation factory function.
  [[nodiscard]] static std::unique_ptr<VirtualAnalysis> factory() { return std::make_unique<SequenceAnalysis>(); }

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

  constexpr static const char PF3D7_IDENT_[]{"Pf3D7_64"};
  std::shared_ptr<const GenomeReference> genome_3D7_ptr_;
  std::shared_ptr<const GenomeCollection> all_reference_genomes_ptr_;
  std::shared_ptr<MutateGenes> mutate_genes_ptr_; // Perform transcript level mutations for all genomes.


  std::string ident_work_directory_;
  constexpr static const char *VARIANT_COUNT_EXT_{".csv"};

};



} // namespace


#endif //KGL_KGL_ANALYSIS_CHECK_H