//
// Created by kellerberrin on 4/6/20.
//

#ifndef KGL_ANALYSIS_INFO_FILTER_H
#define KGL_ANALYSIS_INFO_FILTER_H


#include "kgl_runtime.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_unphased_population.h"
#include "kgl_analysis_null.h"
#include "kgl_filter.h"


namespace kellerberrin::genome {   //  organization::project level namespace

// Test the Info data filters.
class InfoFilterAnalysis : public NullAnalysis {

public:

  InfoFilterAnalysis() = default;
  ~InfoFilterAnalysis() = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "INFO_FILTER"; }
  [[nodiscard]] std::unique_ptr<NullAnalysis> factory() const override { return std::make_unique<InfoFilterAnalysis>(); }

  // Setup the analytics to process VCF data.
  // This function must be redefined.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const RuntimeParameterMap& named_parameters,
                                         std::shared_ptr<const GenomeCollection> reference_genomes) override;

  // Perform the genetic analysis per VCF file (need not be redefined)
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const UnphasedPopulation> vcf_population) override;

  // Perform the genetic analysis per iteration (need not be redefined)
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results (need not be redefined)
  [[nodiscard]] bool finalizeAnalysis() override;

private:


};



} // namespace


#endif //KGL_ANALYSIS_INFO_FILTER_H
