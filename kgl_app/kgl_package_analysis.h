//
// Created by kellerberrin on 4/5/20.
//

#ifndef KGL_PACKAGE_ANALYSIS_H
#define KGL_PACKAGE_ANALYSIS_H


#include "kel_exec_env.h"
#include "kgl_runtime_config.h"
#include "kgl_genome_collection.h"
#include "kgl_variant_db_population.h"
#include "kgl_analysis_all.h" // Includes all the defined active analysis objects.


namespace kellerberrin::genome {   //  organization::project level namespace


// Adds flag to show if analysis is active.
// Analytics that encounter error states disable themselves via this flag.
using VirtualAnalysisArray = std::vector<std::pair<std::unique_ptr<VirtualAnalysis>, bool>>;

// Manages runtime analytics within an execution package.
class PackageAnalysis {

public:

  explicit PackageAnalysis(const RuntimeConfiguration& runtime_contig) : runtime_contig_(runtime_contig) {
    // Defined in "kgl_analysis_all.h"
    registered_analysis_ = getAnalysisVector();

  }

  ~PackageAnalysis() = default;

  // Setup the analytics to process data.
  [[nodiscard]] bool initializeAnalysis(const RuntimePackage& package,
                                        std::shared_ptr<const GenomeCollection> reference_genomes) const;

  // Perform the genetic analysis per VCF file read.
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> file_data) const;

  // Perform the genetic analysis per iteration (multiple files grouped together).
  [[nodiscard]] bool iterationAnalysis() const;

  // All data has been presented, finalize analysis and write results.
  [[nodiscard]] bool finalizeAnalysis() const;

private:

  const RuntimeConfiguration runtime_contig_;
  // All available analytics
  VirtualAnalysisVector registered_analysis_;
  // Active analytics for this package
  mutable VirtualAnalysisArray active_analysis_;


};




} // namespace

#endif //KGL_PACKAGE_ANALYSIS_H
