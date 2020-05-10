//
// Created by kellerberrin on 1/5/20.
//

#include "kgl_package.h"
#include "kgl_variant_factory_vcf.h"


namespace kgl = kellerberrin::genome;


void kgl::ExecutePackage::executeAll() const {

  // for all packages.
  for (auto const& [package_ident, package] : package_map_) {

    // Get reference genomes.
    ExecEnv::log().info("ExecutePackage::executeAll, Load Reference Genomes for Package: {}", package_ident);
    std::shared_ptr<GenomeCollection> reference_genome_ptr = loadReferenceGenomes(package);
    ExecEnv::log().info("ExecutePackage::executeAll, Data Files and perform Analysis for Package: {}", package_ident);

    // Setup the analytics
    if (not package_analysis_.initializeAnalysis(package, reference_genome_ptr)) {

      ExecEnv::log().error("ExecutePackage::executeAll, Problem initializing Analysis for Package: {}", package_ident);

    }

    // Iterate through the VCF files and update analytics.
    for (auto const& iterative_files : package.iterativeFileList()) {

      std::shared_ptr<UnphasedPopulation> vcf_iterative_data = iterateVCFDataFiles(package, reference_genome_ptr, iterative_files);

      if (not package_analysis_.iterateAnalysis(reference_genome_ptr, vcf_iterative_data)) {

        ExecEnv::log().error("ExecutePackage::executeAll, Problem performing Analysis for Package: {}", package_ident);

      }

    }

    // Complete and write the analytics.
    if (not package_analysis_.finalizeAnalysis(reference_genome_ptr)) {

      ExecEnv::log().error("ExecutePackage::executeAll, Problem finalizing Analysis for Package: {}", package_ident);

    }

  }

}


void kgl::ExecutePackage::verifyPackages() const {

  // for all packages.
  for (auto const& [package_ident, package] : package_map_) {

    // Confirm that requested analytics are defined.
    for (auto const& analysis_ident : package.analysisList()) {

      auto result = analysis_map_.find(analysis_ident);
      if (result == analysis_map_.end()) {

        ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Analysis: {}, not defined", package_ident, analysis_ident);

      }

    }

    //confirm that all reference genomes exist
    for (auto const& genome_ident : package.genomeDatabaseList()) {

      auto result = genome_map_.find(genome_ident);
      if (result == genome_map_.end()) {

        ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Reference Genome: {}, not defined", package_ident, genome_ident);

      }

    }


    //confirm that all iterative load files exist.
    // Note that iterativeFileList() returns a nested vector, std::vector<std::vector<std::string>>
    for (auto const& vcf_file_vector : package.iterativeFileList()) {

      for (auto const& vcf_file_ident : vcf_file_vector) {

        auto result = vcf_file_map_.find(vcf_file_ident);
        if (result == vcf_file_map_.end()) {

          ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Iterative load file: {}, not defined", package_ident, vcf_file_ident);

        }

      }

    }

    ExecEnv::log().info("ExecutePackage::verifyPackage, Package: {}, All Reference Genomes, data files and analysis types are defined.", package_ident);

  }

}

std::unique_ptr<kgl::GenomeCollection> kgl::ExecutePackage::loadReferenceGenomes(const RuntimePackage& package) const {

  std::unique_ptr<GenomeCollection> genome_collection_ptr(std::make_unique<GenomeCollection>());

  for (auto const& genome_ident  :  package.genomeDatabaseList()) {

    auto result = genome_map_.find(genome_ident);
    if (result == genome_map_.end()) {

      ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Reference Genome: {}, not defined", package.packageIdentifier(), genome_ident);

    }

     // Create the genome database.
     std::shared_ptr<GenomeReference> genome_ptr = kgl::GenomeReference::createGenomeDatabase(result->second.genomeIdentifier(),
                                                                                               result->second.fastaFileName(),
                                                                                               result->second.gffFileName(),
                                                                                               result->second.gafFileName(),
                                                                                               result->second.translationTable());

    if (not genome_collection_ptr->addGenome(genome_ptr)) {

      ExecEnv::log().error("ExecutePackage::loadReferenceGenomes; Unable to add Genome Database: {} (probable duplicate)", genome_ptr->genomeId());

    }

  }

  return genome_collection_ptr;

}



std::unique_ptr<kgl::UnphasedPopulation> kgl::ExecutePackage::iterateVCFDataFiles( const RuntimePackage& package,
                                                                                   std::shared_ptr<const GenomeCollection> reference_genomes,
                                                                                   const std::vector<std::string>& iterative_files) const {

  std::unique_ptr<UnphasedPopulation> population_ptr(std::make_unique<UnphasedPopulation>(package.packageIdentifier()));

  ExecEnv::log().info("ExecutePackage::loadVCFDataFiles, Package: {}, Iterative File Vector Size:{}", package.packageIdentifier(), iterative_files.size());

  for (auto const& loadfile : iterative_files) {

    ExecEnv::log().info("ExecutePackage::loadVCFDataFiles, Package: {}, VCF file ident: {}", package.packageIdentifier(), loadfile);

    auto result = vcf_file_map_.find(loadfile);
    if (result == vcf_file_map_.end()) {

      ExecEnv::log().critical("ExecutePackage::loadVCFDataFiles, Package: {}, VCF file ident: {}, not defined", package.packageIdentifier(), loadfile);

    }

    std::optional<std::shared_ptr<const GenomeReference>> ref_genome_opt = reference_genomes->getOptionalGenome(result->second.referenceGenome());

    if (not ref_genome_opt) {

      ExecEnv::log().critical("Package: {}, Reference Genome {} Not Found for VCF file ident: {}",
                              package.packageIdentifier(), result->second.referenceGenome(), loadfile);

    }

    if (result->second.parserType() == VCFParserEnum::GRChNoGenome) {


      // Read variants.
      std::shared_ptr<UnphasedGenome> parsed_variants = VcfFactory::GRChNoGenomeVCFVariants(ref_genome_opt.value(),
                                                                                            result->second.fileName(),
                                                                                            contig_alias_);

      std::pair<size_t, size_t> valid_count = parsed_variants->validate(ref_genome_opt.value());
      ExecEnv::log().info("Genome: {}, Total Variants: {}, Validated Variants: {}", parsed_variants->genomeId(), valid_count.first, valid_count.second);

      std::pair<size_t, size_t> merge_stats = population_ptr->mergeUniqueGenome(parsed_variants);
      ExecEnv::log().info("Population: {} merges Genome: {}, Variants Presented: {}, Unique Variants Accepted: {}",
                           population_ptr->populationId(), parsed_variants->genomeId(), merge_stats.first, merge_stats.second);

    } else if (result->second.parserType() == VCFParserEnum::GatkMultiGenome) {
// Pfalciparum VCF files handled here.

      // Read variants.
      std::shared_ptr<UnphasedPopulation> parsed_variants = VcfFactory::gatkMultiGenomeVCFVariants( ref_genome_opt.value(), result->second.fileName());

      std::pair<size_t, size_t> valid_count = parsed_variants->validate(ref_genome_opt.value());
      ExecEnv::log().info("Population: {}, Total Variants: {}, Validated Variants: {}", parsed_variants->populationId(), valid_count.first, valid_count.second);

      population_ptr->mergePopulation(parsed_variants);
      ExecEnv::log().info("Population: {} merges {} Variants, total Variants: {}", population_ptr->populationId(), parsed_variants->variantCount(), population_ptr->variantCount());

    } else {

      ExecEnv::log().critical("ExecutePackage::loadVCFDataFiles, Package: {}, VCF file ident: {}, parser not implemented (only Gatk, GRCh available)", package.packageIdentifier(), loadfile);

    }

  }

  return population_ptr;

}

