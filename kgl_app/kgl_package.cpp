//
// Created by kellerberrin on 1/5/20.
//

#include "kgl_package.h"
#include "kgl_variant_factory_parsers.h"
#include "kol_OntologyDatabase.h"

namespace kgl = kellerberrin::genome;
namespace kol = kellerberrin::ontology;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kgl::ExecutePackage::executeActive() const {

  // for all active packages.
  for (auto const& active_package : runtime_config_.activePackages()) {

    // Get package definition
    auto result = runtime_config_.runtimePackageMap().find(active_package.packageIdentifier());

    if (result == runtime_config_.runtimePackageMap().end()) {

      ExecEnv::log().error( "ExecutePackage::executeActive, Active package :{} could not find matching package definition",
                            active_package.packageIdentifier());

      continue; // execute next active package.

    }

    auto [package_ident, package] = *result;

    // Get reference genomes.
    ExecEnv::log().info("Load Runtime Resources for Package: {}", package_ident);
    std::shared_ptr<GenomeCollection> reference_genome_ptr = loadRuntimeResources(package);
    ExecEnv::log().info("Load Data Files and perform Analysis for Package: {}", package_ident);

    // Setup the analytics
    if (not package_analysis_.initializeAnalysis(package, reference_genome_ptr)) {

      ExecEnv::log().error("ExecutePackage::executeActive, Problem initializing Analysis for Package: {}", package_ident);

    }

    // Iterate through the VCF files and update analytics.
    for (auto const& iterative_files : package.iterativeFileList()) {

      for (auto const& data_file : iterative_files) {

        std::shared_ptr<DataDB> data_ptr = readDataFile(package, reference_genome_ptr, data_file);

        if (not package_analysis_.fileReadAnalysis(data_ptr)) {

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


std::unique_ptr<kgl::GenomeCollection> kgl::ExecutePackage::loadRuntimeResources(const RuntimePackage& package) const {

  std::unique_ptr<GenomeCollection> genome_collection_ptr(std::make_unique<GenomeCollection>());

  for (auto const& [resource_ident, resource_type]  :  package.resourceDatabaseList()) {

    if (resource_type == RuntimeResourceType::GENOME_DATABASE) {


      auto result = runtime_config_.resourceMap().find(resource_ident);
      if (result == runtime_config_.resourceMap().end()) {

        ExecEnv::log().critical("ExecutePackage::loadRuntimeResources, Package: {}, Reference Genome: {}, not defined", package.packageIdentifier(), resource_ident);

      }

      auto const& [resource_ident, resource_base_ptr] = *result;
      auto genome_resource_ptr = std::dynamic_pointer_cast<const RuntimeGenomeResource>(resource_base_ptr);

      if (not genome_resource_ptr) {

        ExecEnv::log().critical("ExecutePackage::loadRuntimeResources, Resource: {} is not a Genome Database", resource_ident);

      }

      // Create the genome database.
      std::shared_ptr<GenomeReference> genome_ptr = kgl::GenomeReference::createGenomeDatabase(genome_resource_ptr->genomeIdentifier(),
                                                                                               genome_resource_ptr->fastaFileName(),
                                                                                               genome_resource_ptr->gffFileName(),
                                                                                               genome_resource_ptr->gafFileName(),
                                                                                               genome_resource_ptr->idFileName(),
                                                                                               genome_resource_ptr->translationTable());

      if (not genome_collection_ptr->addGenome(genome_ptr)) {

        ExecEnv::log().error("ExecutePackage::loadRuntimeResources; Unable to add Genome Database: {} (probable duplicate)", genome_ptr->genomeId());

      }

    } else if (resource_type == RuntimeResourceType::ONTOLOGY_DATABASE) {

      auto result = runtime_config_.resourceMap().find(resource_ident);
      if (result == runtime_config_.resourceMap().end()) {

        ExecEnv::log().critical("ExecutePackage::loadRuntimeResources, Package: {}, Ontology Database: {}, not defined", package.packageIdentifier(), resource_ident);

      }

      auto const& [resource_ident, resource_base_ptr] = *result;
      auto ontology_resource_ptr = std::dynamic_pointer_cast<const RuntimeOntologyResource>(resource_base_ptr);

      if (not ontology_resource_ptr) {

        ExecEnv::log().critical("ExecutePackage::loadRuntimeResources, Resource: {} is not an Ontology Database", resource_ident);

      }

      std::shared_ptr<const kol::OntologyDatabase> ontology_ptr(std::make_shared<const kol::OntologyDatabase>( ontology_resource_ptr->goGraphFileName(),
                                                                                                               ontology_resource_ptr->annotationFileName()));

      ExecEnv::log().info("ExecutePackage::loadRuntimeResources; **********Creating Ontology Database: {}", resource_ident);
      ExecEnv::log().info("ExecutePackage::loadRuntimeResources; **********Gaf File : {}, go terms : {}, Genes: {}",
                          ontology_resource_ptr->annotationFileName(), ontology_ptr->annotation()->getNumGoTerms(), ontology_ptr->annotation()->getNumGenes());
      ExecEnv::log().info("ExecutePackage::loadRuntimeResources; ***********Go File: {}, Vertices: {}, Edges: {}",
                          ontology_resource_ptr->goGraphFileName(), ontology_ptr->goGraph()->getNumVertices(), ontology_ptr->goGraph()->getNumEdges());

    }

  }

  return genome_collection_ptr;

}


std::shared_ptr<kgl::DataDB>
kgl::ExecutePackage::readDataFile(const RuntimePackage& package,
                                  std::shared_ptr<const GenomeCollection> reference_genomes,
                                  const std::string& data_file) const {


  ExecEnv::log().info("Package: {}, Data file ident: {}", package.packageIdentifier(), data_file);

  auto result = runtime_config_.dataFileMap().find(data_file);
  if (result == runtime_config_.dataFileMap().end()) {

    ExecEnv::log().critical("ExecutePackage::readDataFile, Package: {}, data file ident: {}, not defined",
                             package.packageIdentifier(), data_file);

  }

  auto [file_ident, file_info_ptr] = *result;

  // Selects the appropriate parser and returns a base class data object.
  std::shared_ptr<kgl::DataDB> data_ptr = ParserSelection::parseData(reference_genomes, file_info_ptr, runtime_config_.evidenceMap(), runtime_config_.contigAlias());

  return data_ptr;

}

