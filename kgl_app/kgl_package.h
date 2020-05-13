//
// Created by kellerberrin on 1/5/20.
//

#ifndef KGL_PACKAGE_H
#define KGL_PACKAGE_H

#include "kel_exec_env.h"
#include "kgl_runtime.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_unphased_population.h"
#include "kgl_package_analysis.h"

namespace kellerberrin::genome {   //  organization::project level namespace

class ExecutePackage {

public:


  ExecutePackage( const RuntimeProperties& runtime_options, const std::string& work_directory)
                  : contig_alias_(runtime_options.getContigAlias()),
                    vcf_file_map_(runtime_options.getVCFFiles()),
                    genome_map_(runtime_options.getGenomeReferenceMap()),
                    analysis_map_(runtime_options.getAnalysisMap()),
                    package_map_(runtime_options.getPackageMap()),
                    evidence_map_(runtime_options.getEvidenceMap()),
                    package_analysis_(work_directory, analysis_map_) { verifyPackages(); }

  ~ExecutePackage() = default;

  // Executes all the application logic.
  void executeAll() const;

private:

  // The Runtime information loaded from the XML file.
  const ContigAliasMap contig_alias_;
  const RuntimeVCFFileMap vcf_file_map_;
  const RuntimeGenomeDatabaseMap genome_map_;
  const RuntimeAnalysisMap analysis_map_;
  const RuntimePackageMap package_map_;
  const VariantEvidenceMap evidence_map_;
  // The analysis management object.
  const PackageAnalysis package_analysis_;

  // Check the integrity of all the XML information.
  void verifyPackages() const;
  // Load the reference genomes.
  [[nodiscard]] std::unique_ptr<GenomeCollection> loadReferenceGenomes(const RuntimePackage& package) const;
  // Load the VCF data per iteration
  [[nodiscard]] std::unique_ptr<UnphasedPopulation> iterateVCFDataFiles( const RuntimePackage& package,
                                                                         std::shared_ptr<const GenomeCollection> reference_genomes,
                                                                         const std::vector<std::string>& iterative_files) const;
};


} // namespace

#endif //KGL_PACKAGE_H
