//
// Created by kellerberrin on 1/5/20.
//

#include "kgl_package.h"
#include "kgl_ped_parser.h"
#include "kgl_variant_factory_pf3k_impl.h"
#include "kgl_variant_factory_grch_impl.h"
#include "kgl_variant_factory_1000_impl.h"


namespace kgl = kellerberrin::genome;


void kgl::ExecutePackage::executeActive() const {

  // for all active packages.
  for (auto const& active_package : active_packages_) {

    // Get package definition
    auto result = package_map_.find(active_package.packageIdentifier());

    if (result == package_map_.end()) {

      ExecEnv::log().error( "ExecutePackage::executeActive, Active package :{} could not find matching package definition",
                            active_package.packageIdentifier());

      continue; // execute next active package.

    }

    auto [package_ident, package] = *result;

    // Get reference genomes.
    ExecEnv::log().info("Load Reference Genomes for Package: {}", package_ident);
    std::shared_ptr<GenomeCollection> reference_genome_ptr = loadReferenceGenomes(package);
    ExecEnv::log().info("Load Data Files and perform Analysis for Package: {}", package_ident);

    // Setup the analytics
    if (not package_analysis_.initializeAnalysis(package, reference_genome_ptr)) {

      ExecEnv::log().error("ExecutePackage::executeActive, Problem initializing Analysis for Package: {}", package_ident);

    }

    // Iterate through the VCF files and update analytics.
    for (auto const& iterative_files : package.iterativeFileList()) {

      for (auto const& data_file : iterative_files) {

        std::shared_ptr<DataObjectBase> read_data_object = readDataFiles(package, reference_genome_ptr, data_file);

        if (not package_analysis_.fileReadAnalysis(read_data_object)) {

          ExecEnv::log().error("ExecutePackage::executeActive, Problem performing Read File Analysis for Package: {}", package_ident);

        }

      }

      if (not package_analysis_.iterationAnalysis()) {

        ExecEnv::log().error("ExecutePackage::executeActive, Problem performing Analysis for Package: {}", package_ident);

      }

    }

    // Complete and write the analytics.
    if (not package_analysis_.finalizeAnalysis()) {

      ExecEnv::log().error("ExecutePackage::executeActive, Problem finalizing Analysis for Package: {}", package_ident);

    }

  }

}


void kgl::ExecutePackage::verifyPackages() const {

  // List the active packages.
  for (auto const& active : active_packages_) {

    ExecEnv::log().info("Sequentially Executing Active Package: {}", active.packageIdentifier());

  }

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

        auto result = data_file_map_.find(vcf_file_ident);
        if (result == data_file_map_.end()) {

          ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Iterative load file: {}, not defined", package_ident, vcf_file_ident);

        }

      }

    }

    ExecEnv::log().info("Package: {}, All Reference Genomes, data files and analysis types are defined.", package_ident);

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


std::shared_ptr<kgl::DataObjectBase>
kgl::ExecutePackage::readDataFiles(const RuntimePackage& package,
                                   std::shared_ptr<const GenomeCollection> reference_genomes,
                                   const std::string& data_file) const {


  ExecEnv::log().info("Package: {}, Data file ident: {}", package.packageIdentifier(), data_file);

  auto result = data_file_map_.find(data_file);
  if (result == data_file_map_.end()) {

    ExecEnv::log().critical("ExecutePackage::readDataFiles, Package: {}, data file ident: {}, not defined",
                             package.packageIdentifier(), data_file);

  }

  auto [file_ident, file_info_ptr] = *result;

  switch(file_info_ptr->parserType()) {

    case DataFileParserEnum::GatkMultiGenome:
      return readVCF<Pf3kVCFImpl, UnphasedPopulation>(reference_genomes, file_info_ptr);

    case DataFileParserEnum::GRChNoGenome:
      return readVCF<GrchVCFImpl, UnphasedPopulation>(reference_genomes, file_info_ptr);

    case DataFileParserEnum::MultiGenomePhased:
      return readVCF<Genome1000VCFImpl, DiploidPopulation>(reference_genomes, file_info_ptr);

    case DataFileParserEnum::PedAncestry:
      return readPEDAncestry(file_info_ptr);

    case DataFileParserEnum::NotImplemented:
      ExecEnv::log().critical("ExecutePackage::readDataFiles, Package: {}, data file ident: {}, parser not defined",
                               package.packageIdentifier(), data_file);

  }

  return readPEDAncestry(file_info_ptr); // Never reached.

}


// Read and parse the specified ancestry file.
std::shared_ptr<kgl::DataObjectBase> kgl::ExecutePackage::readPEDAncestry(std::shared_ptr<BaseFileInfo> file_info) const {

  auto ped_file_info = std::dynamic_pointer_cast<PedAncestryInfo>(file_info);

  if (not ped_file_info) {

    ExecEnv::log().critical("ExecutePackage::readDataFiles, Expected PED (ancestry .ped) file for file ident: {}", file_info->identifier());

  }

  std::shared_ptr<GenomePEDData> ped_data(std::make_shared<GenomePEDData>(ped_file_info->identifier()));

  ParsePedFile ped_parser(ped_data);

  ped_parser.readParsePEDImpl(ped_file_info->fileName());

  return ped_data;

}