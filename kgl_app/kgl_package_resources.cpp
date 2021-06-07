//
// Created by kellerberrin on 14/4/21.
//

#include "kgl_package.h"
#include "kgl_variant_factory_parsers.h"
#include "kgl_ontology_database.h"

namespace kgl = kellerberrin::genome;



std::shared_ptr<const kgl::AnalysisResources> kgl::ExecutePackage::loadRuntimeResources(const RuntimePackage& package) const {

  std::shared_ptr<AnalysisResources> resource_ptr(std::make_shared<AnalysisResources>());

  for (auto const& [resource_type, resource_ident]  :  package.resourceDatabaseList()) {

    switch (resource_type) {

      case RuntimeResourceType::GENOME_DATABASE:
        loadGenomeResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::ONTOLOGY_DATABASE:
        loadOntologyResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::GENE_NOMENCLATURE:
        loadGeneNomenclatureResource(resource_ident, resource_ptr);
        break;

      default:
        ExecEnv::log().critical("ExecutePackage::loadRuntimeResources, Package: {} Attempt to load unknown resource type, resource ident: {}.",
                                package.packageIdentifier(), resource_ident);
        break;

    }

  }

  return resource_ptr;

}


void kgl::ExecutePackage::loadGenomeResource(const std::string& genome_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {


  auto result = runtime_config_.resourceMap().find(genome_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadGenomeResource, Reference Genome: {}, not defined", genome_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto genome_resource_ptr = std::dynamic_pointer_cast<const RuntimeGenomeResource>(resource_base_ptr);

  if (not genome_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadGenomeResource, Resource: {} is not a Genome Database", resource_ident);

  }

  // Create the genome database.
  std::shared_ptr<GenomeReference> genome_ptr = kgl::GenomeReference::createGenomeDatabase(genome_resource_ptr->genomeIdentifier(),
                                                                                           genome_resource_ptr->fastaFileName(),
                                                                                           genome_resource_ptr->gffFileName(),
                                                                                           genome_resource_ptr->gafFileName(),
                                                                                           genome_resource_ptr->idFileName(),
                                                                                           genome_resource_ptr->translationTable());

  resource_ptr->addResource(genome_ptr);

}

void kgl::ExecutePackage::loadOntologyResource(const std::string& ontology_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(ontology_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadOntologyResource, Ontology Database: {}, not defined", ontology_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto ontology_resource_ptr = std::dynamic_pointer_cast<const RuntimeOntologyResource>(resource_base_ptr);

  if (not ontology_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadOntologyResource, Resource: {} is not an Ontology Database", resource_ident);

  }

  std::shared_ptr<const kol::OntologyDatabase> ontology_ptr(std::make_shared<const kol::OntologyDatabase>( ontology_ident,
                                                                                                           ontology_resource_ptr->goGraphFileName(),
                                                                                                           ontology_resource_ptr->annotationFileName()));

  resource_ptr->addResource(ontology_ptr);

}

void kgl::ExecutePackage::loadGeneNomenclatureResource(const std::string& nomenclature_ident, const std::shared_ptr<AnalysisResources>&) const {

  ExecEnv::log().info("ExecutePackage::loadGeneNomenclatureResource; **********Loading Gene NOMENCLATURE resource: {}", nomenclature_ident);

}
