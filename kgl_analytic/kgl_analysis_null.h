//
// Created by kellerberrin on 4/5/20.
//

#ifndef KGL_NULL_ANALYSIS_H
#define KGL_NULL_ANALYSIS_H


#include "kgl_runtime.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_unphased_population.h"



namespace kellerberrin::genome {   //  organization::project level namespace

// The base class for all analytics.
class NullAnalysis {

public:

  NullAnalysis() = default;
  ~NullAnalysis() = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] virtual std::string ident() const { return "NULL"; }
  [[nodiscard]] virtual std::unique_ptr<NullAnalysis> factory() const { return std::make_unique<NullAnalysis>(); }

  // Setup the analytics to process VCF data.
  // This function must be redefined.
  [[nodiscard]] virtual bool initializeAnalysis( const std::string& work_directory,
                                                 const RuntimeParameterMap& named_parameters,
                                                 std::shared_ptr<const GenomeCollection> reference_genomes);

  // Perform the genetic analysis per VCF file (need not be redefined)
  [[nodiscard]] virtual bool fileReadAnalysis(std::shared_ptr<const UnphasedPopulation> vcf_iterative_dat);

  // Perform the genetic analysis per iteration (need not be redefined)
  [[nodiscard]] virtual bool iterationAnalysis();

  // All VCF data has been presented, finalize analysis and write results (need not be redefined)
  [[nodiscard]] virtual bool finalizeAnalysis();

private:

};





} // namespace

#endif //KGL_NULL_ANALYSIS_H
