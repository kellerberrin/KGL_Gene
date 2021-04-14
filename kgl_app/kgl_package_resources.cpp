//
// Created by kellerberrin on 14/4/21.
//

#include "kgl_package.h"
#include "kgl_variant_factory_parsers.h"
#include "kol_OntologyDatabase.h"

namespace kgl = kellerberrin::genome;
namespace kol = kellerberrin::ontology;



std::shared_ptr<const kgl::AnalysisResources> kgl::ExecutePackage::loadRuntimeResources(const RuntimePackage& package) const {

  std::shared_ptr<AnalysisResources> resource_ptr(std::make_unique<AnalysisResources>());

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

      if (not resource_ptr->addGenome(genome_ptr)) {

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

    } else if (resource_type == RuntimeResourceType::GENE_NOMENCLATURE) {

      ExecEnv::log().info("ExecutePackage::loadRuntimeResources; **********Loading Gene NOMENCLATURE resource: {}", resource_ident);

    }

  }

  return resource_ptr;

}

