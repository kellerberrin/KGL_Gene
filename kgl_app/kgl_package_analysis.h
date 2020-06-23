//
// Created by kellerberrin on 4/5/20.
//

#ifndef KGL_PACKAGE_ANALYSIS_H
#define KGL_PACKAGE_ANALYSIS_H


#include "kel_exec_env.h"
#include "kgl_runtime.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_unphased_population.h"
#include "kgl_analysis_all.h" // Includes all the defined active analysis objects.


namespace kellerberrin::genome {   //  organization::project level namespace


// Adds flag to show if analysis is active.
// Analytics that encounter error states disable themselves via this flag.
using VirtualAnalysisArray = std::vector<std::pair<std::unique_ptr<VirtualAnalysis>, bool>>;

// Manages runtime analytics within an execution package.
class PackageAnalysis {

public:

  explicit PackageAnalysis(std::string work_directory, const RuntimeAnalysisMap& analysis_map)
  : work_directory_(std::move(work_directory)), analysis_map_(analysis_map) {

    // Defined in "kgl_analysis_all.h"
    registered_analysis_ = getAnalysisVector();

  }
  ~PackageAnalysis() = default;

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis(const RuntimePackage& package,
                                        std::shared_ptr<const GenomeCollection> reference_genomes) const;

  // Perform the genetic analysis per VCF file read.
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const UnphasedPopulation> vcf_iterative_dat) const;

  // Perform the genetic analysis per iteration (multiple VCF files grouped together).
  [[nodiscard]] bool iterationAnalysis() const;

  // All VCF data has been presented, finalize analysis and write results.
  [[nodiscard]] bool finalizeAnalysis() const;

private:

  // Output directory for analysis files.
  const std::string work_directory_;
  // Analysis parameters and details.
  const RuntimeAnalysisMap analysis_map_;
  // All available analytics
  VirtualAnalysisVector registered_analysis_;
  // Active analytics for this package
  mutable VirtualAnalysisArray active_analysis_;


};




} // namespace

#endif //KGL_PACKAGE_ANALYSIS_H
