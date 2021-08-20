//
// Created by kellerberrin on 14/4/21.
//

#include "kgl_package.h"
#include "kgl_variant_factory_parsers.h"
#include "kgl_ontology_database.h"
#include "kgl_uniprot_parser.h"
#include "kgl_citation_parser.h"
#include "kgl_entrez_parser.h"
#include "kgl_bio_pmid_parser.h"



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

      case RuntimeResourceType::ALLELE_CITATION:
        loadCitationResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::ENTREZ_GENE:
        loadEntrezGeneResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::BIO_PMID:
        loadPMIDBioResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::PUBMED_API:
        loadPubmedAPIResource(resource_ident, resource_ptr);
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


void kgl::ExecutePackage::loadCitationResource(const std::string& citation_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(citation_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadCitationResource, Allele Citation Database: {}, not defined", citation_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto citation_resource_ptr = std::dynamic_pointer_cast<const RuntimeCitationResource>(resource_base_ptr);

  if (not citation_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadCitationResource, Resource: {} is not an Allele Citation Resource", resource_ident);

  }

  ParseCitations citation_parser;
  if (not citation_parser.parseCitationFile(citation_resource_ptr->citationFileName())) {

    ExecEnv::log().critical("ExecutePackage::loadCitationResource; failed to create citation resource from file: {}", citation_resource_ptr->citationFileName());

  }

  std::shared_ptr<CitationResource> citation_ptr(std::make_shared<CitationResource>(citation_resource_ptr->citationIdentifier(), citation_parser.getCitationMap()));

  resource_ptr->addResource(citation_ptr);

}


void kgl::ExecutePackage::loadEntrezGeneResource(const std::string& entrez_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(entrez_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadEntrezGeneResource, Entrez Gene Database: {}, not defined", entrez_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto entrez_resource_ptr = std::dynamic_pointer_cast<const RuntimeEntrezResource>(resource_base_ptr);

  if (not entrez_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadEntrezGeneResource, Resource: {} is not an Entrez Gene Resource", resource_ident);

  }

  ParseEntrez entrez_parser;
  if (not entrez_parser.parseEntrezFile(entrez_resource_ptr->entrezFileName())) {

    ExecEnv::log().critical("ExecutePackage::loadEntrezGeneResource; failed to create Entrez resource from file: {}", entrez_resource_ptr->entrezFileName());

  }

  std::shared_ptr<EntrezResource> entrez_ptr(std::make_shared<EntrezResource>(entrez_resource_ptr->entrezIdentifier(), entrez_parser.getEntrezVector()));

  resource_ptr->addResource(entrez_ptr);

}


void kgl::ExecutePackage::loadPMIDBioResource(const std::string& bio_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(bio_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadPMIDBioResource, PMID Bio Database: {}, not defined", bio_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto bio_resource_ptr = std::dynamic_pointer_cast<const RuntimeBioPMIDResource>(resource_base_ptr);

  if (not bio_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadPMIDBioResource, Resource: {} is not an PMID Bio Resource", resource_ident);

  }

  ParseBioPMID bio_pmid_parser;
  if (not bio_pmid_parser.parseBioPMIDRecords(bio_resource_ptr->bioFileName())) {

    ExecEnv::log().critical("ExecutePackage::loadPMIDBioResource; failed to create Bio PMID resource from file: {}", bio_resource_ptr->bioFileName());

  }

  std::shared_ptr<BioPMIDResource> bio_pmid_ptr(std::make_shared<BioPMIDResource>(bio_resource_ptr->bioIdentifier(),
                                                                                  bio_pmid_parser.moveDiseasePMIDMap(),
                                                                                  bio_pmid_parser.moveEntrezPMIDMap()));

  resource_ptr->addResource(bio_pmid_ptr);

}


void kgl::ExecutePackage::loadPubmedAPIResource(const std::string& api_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(api_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadPubmedAPIResource, Pubmed API Database: {}, not defined", api_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto api_resource_ptr = std::dynamic_pointer_cast<const RuntimePubmedAPIResource>(resource_base_ptr);

  if (not api_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadPubmedAPIResource, Resource: {} is not a Pubmed API Resource", resource_ident);

  }

  std::shared_ptr<PubmedRequester> pubmed_api_ptr(std::make_shared<PubmedRequester>(api_resource_ptr->apiIdentifier(), api_resource_ptr->cacheFileName()));

  resource_ptr->addResource(pubmed_api_ptr);

}