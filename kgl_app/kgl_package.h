//
// Created by kellerberrin on 1/5/20.
//

#ifndef KGL_PACKAGE_H
#define KGL_PACKAGE_H

#include "kel_exec_env.h"
#include "kgl_runtime.h"
#include "kgl_genome_db.h"
#include "kgl_package_analysis.h"

namespace kellerberrin::genome {   //  organization::project level namespace

using PopulationPair = std::pair<std::unique_ptr<UnphasedPopulation>, std::unique_ptr<PhasedPopulation>>;

class ExecutePackage {

public:


  ExecutePackage( const RuntimeProperties& runtime_options, const std::string& work_directory)
                  : active_packages_(runtime_options.getActivePackages()),
                    contig_alias_(runtime_options.getContigAlias()),
                    data_file_map_(runtime_options.getDataFiles()),
                    genome_map_(runtime_options.getGenomeReferenceMap()),
                    analysis_map_(runtime_options.getAnalysisMap()),
                    package_map_(runtime_options.getPackageMap()),
                    evidence_map_(runtime_options.getEvidenceMap()),
                    package_analysis_(work_directory, analysis_map_) { verifyPackages(); }

  ~ExecutePackage() = default;

  // Executes all the application logic.
  void executeActive() const;

private:

  // The Runtime information loaded from the XML file.
  const ActivePackageVector active_packages_;
  const ContigAliasMap contig_alias_;
  const RuntimeDataFileMap data_file_map_;
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

  // Load a specified Data file and return a population or ancestry.
  [[nodiscard]] std::shared_ptr<DataObjectBase> readDataFiles(const RuntimePackage& package,
                                                              std::shared_ptr<const GenomeCollection> reference_genomes,
                                                              const std::string& data_file) const;

  // Called when the VCF file specifies the GRChNoGenome parser (gnomad, clinvar, SNPdb etc.)
  [[nodiscard]] std::shared_ptr<DataObjectBase> readGRChNoGenomeVCF(std::shared_ptr<const GenomeCollection> reference_genomes,
                                                                    std::shared_ptr<BaseFileInfo> file_info) const;

  // Called when the VCF file specifies the GatkMultiGenome parser (GatK generated P.Falciparum VCF files)
  [[nodiscard]] std::shared_ptr<DataObjectBase> readGatkMultiGenome(std::shared_ptr<const GenomeCollection> reference_genomes,
                                                                    std::shared_ptr<BaseFileInfo> file_info) const;

  // Called when the VCF file specifies the MultiGenomePhased parser (The 1000 Genomes project)
  [[nodiscard]] std::shared_ptr<DataObjectBase> readMultiGenomePhased(std::shared_ptr<const GenomeCollection> reference_genomes,
                                                                      std::shared_ptr<BaseFileInfo> file_info) const;

  // Called when the PedAncestry parser is specified to read in an ancestry (.ped) file
  [[nodiscard]] std::shared_ptr<DataObjectBase> readPEDAncestry(std::shared_ptr<BaseFileInfo> file_info) const;

};


} // namespace

#endif //KGL_PACKAGE_H
