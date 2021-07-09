//
// Created by kellerberrin on 14/4/21.
//

#include "kgl_package.h"
#include "kgl_variant_factory_parsers.h"
#include "kgl_ontology_database.h"
#include "kgl_uniprot_parser.h"

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
        loadHsGeneNomenclatureResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::GENOME_GENEALOGY:
        loadHsGenomeGenealogyResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::GENOME_AUX_INFO:
        loadHsGenomeAuxResource(resource_ident, resource_ptr);
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

void kgl::ExecutePackage::loadHsGeneNomenclatureResource(const std::string& nomenclature_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(nomenclature_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Nomenclature Database: {}, not defined", nomenclature_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto nomenclature_resource_ptr = std::dynamic_pointer_cast<const RuntimeNomenclatureResource>(resource_base_ptr);

  if (not nomenclature_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadNomenclatureResource, Resource: {} is not a Nomenclature Database", resource_ident);

  }

  if (resource_ident == ResourceBase::NOMENCLATURE_UNIPROTID) {

    ParseUniprotId parse_uniprot;
    if (not parse_uniprot.parseUniprotFile(nomenclature_resource_ptr->nomenclatureFileName())) {

      ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Unable to parse Uniprot Nomenclature file: {}", nomenclature_resource_ptr->nomenclatureFileName());

    }

    std::shared_ptr<const UniprotResource> gene_id_resource(std::make_shared<const UniprotResource>(resource_ident, parse_uniprot.getUniproResource()));

    resource_ptr->addResource(gene_id_resource);

  } else if (resource_ident == ResourceBase::NOMENCLATURE_ENSEMBL) {

    ParseGeneIdents parse_gene_idents;
    if (not parse_gene_idents.parseIdentFile(nomenclature_resource_ptr->nomenclatureFileName())) {

      ExecEnv::log().critical("ExecutePackage::loadGeneNomenclatureResource, Unable to parse Ensembl Nomenclature file: {}", nomenclature_resource_ptr->nomenclatureFileName());

    }

    std::shared_ptr<const EnsemblHGNCResource> gene_id_resource(std::make_shared<const EnsemblHGNCResource>(resource_ident, parse_gene_idents.getSynonymVector()));

    resource_ptr->addResource(gene_id_resource);

  } else {

    ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Nomenclature resource ident: {} unknown", resource_ident);

  }

}


void kgl::ExecutePackage::loadHsGenomeGenealogyResource(const std::string& genealogy_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(genealogy_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadHsGenomeGenealogyResource, Genome Genealogy Database: {}, not defined", genealogy_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto genealogy_resource_ptr = std::dynamic_pointer_cast<const RuntimeGenealogyResource>(resource_base_ptr);

  if (not genealogy_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadHsGenomeGenealogyResource, Resource: {} is not a Genome Genealogy", resource_ident);

  }

  std::shared_ptr<HsGenomeGenealogyData> genealogy_data(std::make_shared<HsGenomeGenealogyData>(genealogy_resource_ptr->genealogyIdentifier()));
  ParseHsGenomeGenealogyFile genealogy_parser(genealogy_data);
  genealogy_parser.readParseImpl(genealogy_resource_ptr->genealogyFileName());

  resource_ptr->addResource(genealogy_data);

}


void kgl::ExecutePackage::loadHsGenomeAuxResource(const std::string& genome_aux_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(genome_aux_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadHsGenomeAuxResource, Genome Genealogy Database: {}, not defined", genome_aux_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto genome_aux_resource_ptr = std::dynamic_pointer_cast<const RuntimeGenomeAuxResource>(resource_base_ptr);

  if (not genome_aux_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadHsGenomeAuxResource, Resource: {} is not a HsGenomeAuxResource", resource_ident);

  }

  std::shared_ptr<HsGenomeAuxData> genome_aux_data(std::make_shared<HsGenomeAuxData>(genome_aux_resource_ptr->genomeAuxIdentifier()));
  ParseHsGenomeAuxFile genome_aux_parser(genome_aux_data);
  genome_aux_parser.readParseImpl(genome_aux_resource_ptr->genomeAuxFileName());

  resource_ptr->addResource(genome_aux_data);

}