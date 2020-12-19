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

  // Called to read VCF files
  template<class VCFParser, class VCFPopulation>
  [[nodiscard]] std::shared_ptr<DataObjectBase> readVCF(std::shared_ptr<const GenomeCollection> reference_genomes,
                                                        std::shared_ptr<BaseFileInfo> file_info) const;

  // Called when the PedAncestry parser is specified to read in an ancestry (.ped) file
  [[nodiscard]] std::shared_ptr<DataObjectBase> readPEDAncestry(std::shared_ptr<BaseFileInfo> file_info) const;

};


template<class VCFParser, class VCFPopulation>
std::shared_ptr<DataObjectBase> ExecutePackage::readVCF( std::shared_ptr<const GenomeCollection> reference_genomes,
                                                         std::shared_ptr<BaseFileInfo> file_info) const {

  auto vcf_file_info = std::dynamic_pointer_cast<RuntimeVCFFileInfo>(file_info);

  if (not vcf_file_info) {

    ExecEnv::log().critical("ExecutePackage::readVCF; Expected VCF file for file ident: {}", file_info->identifier());

  }

  std::optional<std::shared_ptr<const GenomeReference>> ref_genome_opt = reference_genomes->getOptionalGenome(vcf_file_info->referenceGenome());

  if (not ref_genome_opt) {

    ExecEnv::log().critical("ExecutePackage::readVCF; Reference Genome {} Not Found for VCF file ident: {}",
                            vcf_file_info->referenceGenome(), vcf_file_info->identifier());

  }

  auto evidence_opt = evidence_map_.lookupEvidence(vcf_file_info->evidenceIdent());

  if (not evidence_opt) {

    ExecEnv::log().critical("ExecutePackage::readVCF; Evidence Ident {} Not Found for VCF file ident: {}",
                            vcf_file_info->evidenceIdent(), vcf_file_info->identifier());

  }

  // Read variants.
  std::shared_ptr<VCFPopulation> vcf_population_ptr(std::make_shared<VCFPopulation>(vcf_file_info->identifier()));

  VCFParser reader(vcf_population_ptr, ref_genome_opt.value(), contig_alias_, evidence_opt.value());
  reader.readParseVCFImpl(vcf_file_info->fileName());

  auto [total_variants, validated_variants] = vcf_population_ptr->validate(ref_genome_opt.value());

  ExecEnv::log().info("Genome: {}, Total Variants: {}, Validated Variants: {}", vcf_population_ptr->populationId(), total_variants, validated_variants);

  return vcf_population_ptr;

}




} // namespace

#endif //KGL_PACKAGE_H
