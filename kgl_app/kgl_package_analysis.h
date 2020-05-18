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


#include <string>



namespace kellerberrin::genome {   //  organization::project level namespace

// NullAnalysis is the base class, of all available (defined) Analytics.
// Analytic classes are called virtually and provided with data (reference and variant).
using AnalysisVector = std::vector<std::unique_ptr<NullAnalysis>>;
// Active flag. Analytics that encounter error states disable themselves via this flag.
using AnalysisArray = std::vector<std::pair<std::unique_ptr<NullAnalysis>, bool>>; // Adds flag to show if analysis is active.

// Manages runtime analytics within an execution package.
class PackageAnalysis {

public:

  explicit PackageAnalysis(std::string work_directory, const RuntimeAnalysisMap& analysis_map)
  : work_directory_(std::move(work_directory)), analysis_map_(analysis_map) {

    // All available analytics registered here.
    registered_analysis_.push_back(std::make_unique<NullAnalysis>());
    registered_analysis_.push_back(std::make_unique<IntervalAnalysis>());

  }
  ~PackageAnalysis() = default;

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis(const RuntimePackage& package,
                                        std::shared_ptr<const GenomeCollection> reference_genomes) const;

  // Perform the genetic analysis per iteration.
  [[nodiscard]] bool iterateAnalysis(std::shared_ptr<const GenomeCollection> reference_genomes,
                                     std::shared_ptr<const UnphasedPopulation> vcf_iterative_dat) const;

  // All VCF data has been presented, finalize analysis and write results.
  [[nodiscard]] bool finalizeAnalysis(std::shared_ptr<const GenomeCollection> reference_genomes) const;

private:

  // Output directory for analysis files.
  const std::string work_directory_;
  // Analysis parameters and details.
  const RuntimeAnalysisMap analysis_map_;
  // All available analytics
  AnalysisVector registered_analysis_;
  // Active analytics for this package
  mutable  AnalysisArray active_analysis_;


};




} // namespace

#endif //KGL_PACKAGE_ANALYSIS_H
