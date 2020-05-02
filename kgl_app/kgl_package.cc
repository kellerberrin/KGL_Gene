//
// Created by kellerberrin on 1/5/20.
//

#include "kgl_package.h"
#include "kgl_variant_factory_vcf.h"


namespace kgl = kellerberrin::genome;


void kgl::ExecutePackage::executeAll() {

  // for all packages.
  for (auto const& [package_ident, package] : package_map_) {

    ExecEnv::log().info("ExecutePackage::executeAll, Load Reference Genomes for Package: {}", package_ident);

    std::shared_ptr<GenomeCollection> reference_genome_ptr = loadReferenceGenomes(package);
    std::shared_ptr<UnphasedPopulation> vcf_load_data = loadVCFDataFiles(package, reference_genome_ptr);
    std::shared_ptr<UnphasedPopulation> vcf_iterative_data = iterateVCFDataFiles(package, reference_genome_ptr);

  }

}



void kgl::ExecutePackage::verifyPackages() const {

  ExecEnv::log().info("VCF Map Size: {}", vcf_file_map_.size());
  ExecEnv::log().info("Genome Map Size: {}", genome_map_.size());
  ExecEnv::log().info("Analysis Map Size: {}", analysis_map_.size());
  ExecEnv::log().info("Package Map Size: {}", package_map_.size());
  ExecEnv::log().info("Alias Size: {}", contig_alias_.getMap().size());

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

    //confirm that all static load files exist
    for (auto const& vcf_static_ident : package.loadFileList()) {

      auto result = vcf_file_map_.find(vcf_static_ident);
      if (result == vcf_file_map_.end()) {

        ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Static Load file: {}, not defined", package_ident, vcf_static_ident);

      }

    }

    //confirm that all iterative load files exist
    for (auto const& vcf_file_ident : package.iterativeFileList()) {

      auto result = vcf_file_map_.find(vcf_file_ident);
      if (result == vcf_file_map_.end()) {

        ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Iterative load file: {}, not defined", package_ident, vcf_file_ident);

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


std::unique_ptr<kgl::UnphasedPopulation> kgl::ExecutePackage::loadVCFDataFiles( const RuntimePackage& package,
                                                                                std::shared_ptr<const GenomeCollection> reference_genomes) const {

  std::unique_ptr<UnphasedPopulation> population_ptr(std::make_unique<UnphasedPopulation>(package.packageIdentifier()));

  for (auto const& loadfile : package.loadFileList()) {

    auto result = vcf_file_map_.find(loadfile);
    if (result == vcf_file_map_.end()) {

      ExecEnv::log().critical("ExecutePackage::loadVCFDataFiles, Package: {}, VCF file ident: {}, not defined", package.packageIdentifier(), loadfile);

    }

    if (result->second.parserType() == VCFParserEnum::GRChNoGenome) {

      std::optional<std::shared_ptr<const GenomeReference>> ref_genome_opt = reference_genomes->getOptionalGenome(result->second.referenceGenome());

      if (not ref_genome_opt) {

        ExecEnv::log().critical("Package: {}, Reference Genome {} Not Found for VCF file ident: {}",
                                 package.packageIdentifier(), result->second.referenceGenome(), loadfile);

      }

      // Read variants.
      std::shared_ptr<UnphasedGenome> parsed_variants = VcfFactory::GRChNoGenomeVCFVariants(ref_genome_opt.value(),
                                                                                            result->second.fileName(),
                                                                                            contig_alias_);

      std::pair<size_t, size_t> valid_count = parsed_variants->validate(ref_genome_opt.value());
      ExecEnv::log().info("Genome: {}, Total Variants: {}, Validated Variants: {}", parsed_variants->genomeId(), valid_count.first, valid_count.second);

    } else if (result->second.parserType() == VCFParserEnum::GatkMultiGenome) {


    } else {

      ExecEnv::log().critical("ExecutePackage::loadVCFDataFiles, Package: {}, VCF file ident: {}, parser not implemented (only Gatk, GRCh available)", package.packageIdentifier(), loadfile);

    }

  }

  return population_ptr;

}


std::unique_ptr<kgl::UnphasedPopulation> kgl::ExecutePackage::iterateVCFDataFiles( const RuntimePackage& package,
                                                                                   std::shared_ptr<const GenomeCollection> reference_genomes) const {

  ExecEnv::log().info("ExecutePackage::iterateVCFDataFiles, Package: {}, Iterative Files: {},", package.packageIdentifier(), package.iterativeFileList().size());

  std::unique_ptr<UnphasedPopulation> population_ptr(std::make_unique<UnphasedPopulation>(package.packageIdentifier()));

  for (auto const& loadfile : package.iterativeFileList()) {

    ExecEnv::log().info("ExecutePackage::loadVCFDataFiles, Package: {}, VCF file ident: {}", package.packageIdentifier(), loadfile);

    auto result = vcf_file_map_.find(loadfile);
    if (result == vcf_file_map_.end()) {

      ExecEnv::log().critical("ExecutePackage::loadVCFDataFiles, Package: {}, VCF file ident: {}, not defined", package.packageIdentifier(), loadfile);

    }

    if (result->second.parserType() == VCFParserEnum::GRChNoGenome) {

      std::optional<std::shared_ptr<const GenomeReference>> ref_genome_opt = reference_genomes->getOptionalGenome(result->second.referenceGenome());

      if (not ref_genome_opt) {

        ExecEnv::log().critical("Package: {}, Reference Genome {} Not Found for VCF file ident: {}",
                                package.packageIdentifier(), result->second.referenceGenome(), loadfile);

      }

      // Read variants.
      std::shared_ptr<UnphasedGenome> parsed_variants = VcfFactory::GRChNoGenomeVCFVariants(ref_genome_opt.value(),
                                                                                            result->second.fileName(),
                                                                                            contig_alias_);

      std::pair<size_t, size_t> valid_count = parsed_variants->validate(ref_genome_opt.value());
      ExecEnv::log().info("Genome: {}, Total Variants: {}, Validated Variants: {}", parsed_variants->genomeId(), valid_count.first, valid_count.second);

    } else if (result->second.parserType() == VCFParserEnum::GatkMultiGenome) {


    } else {

      ExecEnv::log().critical("ExecutePackage::loadVCFDataFiles, Package: {}, VCF file ident: {}, parser not implemented (only Gatk, GRCh available)", package.packageIdentifier(), loadfile);

    }

  }

  return population_ptr;

}

