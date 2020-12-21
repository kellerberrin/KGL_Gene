//
// Created by kellerberrin on 7/11/20.
//

#ifndef KGL_ANALYSIS_CHECK_H
#define KGL_ANALYSIS_CHECK_H


#include "kgl_sequence/kgl_analysis_virtual.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object checks the underlying population structures for correctness.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization::project level namespace


class VerifyAnalysis : public VirtualAnalysis {

public:

  VerifyAnalysis() = default;
  ~VerifyAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "VERIFY"; }
  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<VerifyAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const RuntimeParameterMap& named_parameters,
                                         std::shared_ptr<const GenomeCollection> reference_genomes) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataObjectBase> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Check for duplicate variants at each location in a contig, in a genome, in a population.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VerifyDuplicates {

public:

  VerifyDuplicates() = default;
  ~VerifyDuplicates() { verifyDuplicates(); }


  bool verifyVariant(const std::shared_ptr<const Variant> variant_ptr);
  size_t duplicateCount() const { return duplicate_count_; }

private:

  std::shared_ptr<const Variant> previous_variant_ptr_;
  std::vector<std::shared_ptr<const Variant>> offset_variants_;
  size_t duplicate_count_{0};

  bool verifyDuplicates();

};


} // namespace


#endif //KGL_KGL_ANALYSIS_CHECK_H
