//
// Created by kellerberrin on 7/11/20.
//

#ifndef KGL_ANALYSIS_CHECK_H
#define KGL_ANALYSIS_CHECK_H


#include "kgl_analysis_virtual.h"

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
                                         const ActiveParameterList& named_parameters,
                                         std::shared_ptr<const GenomeCollection> reference_genomes) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

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
  const std::unordered_map<std::string, size_t>& hashMap() const { return hash_map_; }

private:

  std::shared_ptr<const Variant> previous_variant_ptr_;
  std::vector<std::shared_ptr<const Variant>> offset_variants_;
  size_t duplicate_count_{0};
  std::unordered_map<std::string, size_t> hash_map_;

  bool verifyDuplicates();

};


class VerifyHashDuplicates {

public:

  VerifyHashDuplicates() = default;
  ~VerifyHashDuplicates() = default;


  bool verifyVariant(const std::shared_ptr<const Variant> variant_ptr);
  size_t duplicateCount() const { return duplicate_count_; }
  const std::unordered_map<std::string, size_t>& hashMap() const { return hash_map_; }

private:

  std::unordered_map<std::string, size_t> hash_map_;
  size_t duplicate_count_{0};


};



class VerifyHashFilter {

public:

  VerifyHashFilter() = default;
  ~VerifyHashFilter() = default;


  bool storeVariant(const std::shared_ptr<const Variant> variant_ptr);

  bool checkVariant(const std::shared_ptr<const Variant> variant_ptr);

  size_t errors() const { return errors_; }
  size_t checked() const { return checked_; }

private:

  std::unordered_map<std::string, std::shared_ptr<const Variant>> filtered_map_;
  size_t errors_{0};
  size_t checked_{0};

};




} // namespace


#endif //KGL_KGL_ANALYSIS_CHECK_H
