//
// Created by kellerberrin on 21/6/20.
//

#ifndef KGL_ANALYSIS_VIRTUAL_H
#define KGL_ANALYSIS_VIRTUAL_H



#include "kgl_runtime.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_unphased_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace

// Interface for the analysis class
class VirtualAnalysis;
// Analytic classes are called virtually and provided with data (reference and variant).

class VirtualAnalysis {

public:

  VirtualAnalysis() = default;
  virtual ~VirtualAnalysis() = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] virtual std::string ident() const = 0;
  [[nodiscard]] virtual std::unique_ptr<VirtualAnalysis> factory() const = 0;

  // Setup the analytics to process VCF data.
  [[nodiscard]] virtual bool initializeAnalysis( const std::string& work_directory,
                                                 const RuntimeParameterMap& named_parameters,
                                                 std::shared_ptr<const GenomeCollection> reference_genomes) = 0;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] virtual bool fileReadAnalysis(std::shared_ptr<const UnphasedPopulation> vcf_iterative_dat) = 0;

  // Perform the genetic analysis per iteration
  [[nodiscard]] virtual bool iterationAnalysis() = 0;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] virtual bool finalizeAnalysis() = 0;

private:

};


}

#endif //KGL_KGL_ANALYSIS_VIRT_H
