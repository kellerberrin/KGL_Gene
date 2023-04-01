//
// Created by kellerberrin on 14/4/21.
//

#include "kel_utility.h"

#include "kgl_package.h"
#include "kgl_variant_factory_parsers.h"
#include "kgl_ontology_database.h"
#include "kgl_uniprot_parser.h"
#include "kgl_citation_parser.h"
#include "kgl_entrez_parser.h"
#include "kgl_bio_pmid_parser.h"
#include "kgl_pf7_sample_parser.h"


namespace kgl = kellerberrin::genome;



std::shared_ptr<const kgl::AnalysisResources> kgl::ExecutePackage::loadRuntimeResourcesDef(const RuntimePackage& package) const {

  std::shared_ptr<AnalysisResources> resource_ptr(std::make_shared<AnalysisResources>());

  for (auto const& [resource_type, resource_ident]  :  package.resourceDatabaseDef()) {

    switch (Utility::hash(resource_type)) {

      case Utility::hash(RuntimeProperties::GENOME_DATABASE_):
        loadGenomeResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(RuntimeProperties::ONTOLOGY_DATABASE_):
        loadOntologyResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(RuntimeProperties::GENE_ID_DATABASE_):
        loadHsGeneNomenclatureResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(RuntimeProperties::GENEALOGY_ID_DATABASE_):
        loadHsGenomeGenealogyResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(RuntimeProperties::AUX_ID_DATABASE_):
        loadHsGenomeAuxResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(RuntimeProperties::CITATION_DATABASE_):
        loadCitationResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(RuntimeProperties::ENTREZ_DATABASE_):
        loadEntrezGeneResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(RuntimeProperties::PUBMED_LIT_API_):
        loadPubmedAPIResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(RuntimeProperties::PF7_SAMPLE_DATABASE_):
        loadPf7SampleResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(RuntimeProperties::COMMENT_):
      case Utility::hash(RuntimeProperties::HELP_):
        break;

      default:
        ExecEnv::log().critical("ExecutePackage::loadRuntimeResourcesDef, Package: {} Attempt to load unknown resource type: '{}', resource ident: '{}'.",
                                package.packageIdentifier(), resource_type, resource_ident);
        break;

    }

  }

  return resource_ptr;

}


std::shared_ptr<const kgl::AnalysisResources> kgl::ExecutePackage::loadRuntimeResources(const RuntimePackage& package) const {

  std::shared_ptr<AnalysisResources> resource_ptr(std::make_shared<AnalysisResources>());

  for (auto const& [resource_type, resource_ident]  :  package.resourceDatabaseList()) {

    switch (resource_type) {

      case RuntimeResourceType::GENOME_DATABASE:
//        loadGenomeResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::HSAPIEN_ONTOLOGY:
//        loadOntologyResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::GENE_NOMENCLATURE:
//        loadHsGeneNomenclatureResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::GENOME_GENEALOGY:
//        loadHsGenomeGenealogyResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::GENOME_AUX_INFO:
//        loadHsGenomeAuxResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::ALLELE_CITATION:
//        loadCitationResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::ENTREZ_GENE:
//        loadEntrezGeneResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::PUBMED_API:
//        loadPubmedAPIResource(resource_ident, resource_ptr);
        break;

      case RuntimeResourceType::PF7_SAMPLE_DATA:
//        loadPf7SampleResource(resource_ident, resource_ptr);
        break;

      default:
        ExecEnv::log().critical("ExecutePackage::loadRuntimeResources, Package: {} Attempt to load unknown resource type, resource ident: {}.",
                                package.packageIdentifier(), resource_ident);
        break;

    }

  }

  return resource_ptr;

}


void kgl::ExecutePackage::loadGenomeResource(const std::string& genome_type,
                                             const std::string& genome_ident,
                                             const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(genome_type, genome_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadGenomeResource, Reference Genome: {}, not defined", genome_ident);

  }

  auto const& params = params_opt.value();
  auto fasta_opt = params.getParameter(RuntimeProperties::FASTA_FILE_);
  if (not fasta_opt) {

    ExecEnv::log().critical("ExecutePackage::loadGenomeResource, Reference Genome: {}, Fasta file not defined", genome_ident);

  }
  auto gff_opt = params.getParameter(RuntimeProperties::GFF_FILE_);
  if (not gff_opt) {

    ExecEnv::log().critical("ExecutePackage::loadGenomeResource, Reference Genome: {}, Gff file not defined", genome_ident);

  }
  auto translation_opt = params.getParameter(RuntimeProperties::TRANSLATION_TABLE_);
  if (not translation_opt) {

    ExecEnv::log().critical("ExecutePackage::loadGenomeResource, Reference Genome: {}, Translation table not defined", genome_ident);

  }
  // The Gaf file parameter is optional
  auto gaf_opt = params.getParameter(RuntimeProperties::GAF_ANNOTATION_FILE_);
  std::string gaf_file;
  if (gaf_opt) {

    gaf_file = gaf_opt.value();

  }
  // Create the genome database.
  std::shared_ptr<GenomeReference> genome_ptr = kgl::GenomeReference::createGenomeDatabase(genome_ident,
                                                                                           fasta_opt.value(),
                                                                                           gff_opt.value(),
                                                                                           gaf_file,
                                                                                           translation_opt.value());

  resource_ptr->addResource(genome_ptr);

}

void kgl::ExecutePackage::loadOntologyResource(const std::string& resource_type,
                                               const std::string& ontology_ident,
                                               const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, ontology_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadOntologyResource, Ontology Database: {}, not defined", ontology_ident);

  }

  auto const& params = params_opt.value();
  auto go_graph_opt = params.getParameter(RuntimeProperties::GO_ONTOLOGY_FILE_);
  if (not go_graph_opt) {

    ExecEnv::log().critical("ExecutePackage::loadOntologyResource, Ontology Ident: {} Go graph file not defined", ontology_ident);

  }

  auto gaf_opt = params.getParameter(RuntimeProperties::GAF_ANNOTATION_FILE_);
  if (not gaf_opt) {

    ExecEnv::log().critical("ExecutePackage::loadOntologyResource, Ontology Ident: {} Gaf file not defined", ontology_ident);

  }

  std::shared_ptr<const kol::OntologyDatabase> ontology_ptr(std::make_shared<const kol::OntologyDatabase>( ontology_ident,
                                                                                                           go_graph_opt.value(),
                                                                                                           gaf_opt.value()));

  resource_ptr->addResource(ontology_ptr);

}

void kgl::ExecutePackage::loadHsGeneNomenclatureResource(const std::string& resource_type,
                                                         const std::string& nomenclature_ident,
                                                         const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, nomenclature_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Nomenclature Database: {}, not defined", nomenclature_ident);

  }

  auto const& params = params_opt.value();
  auto file_name_opt = params.getParameter(RuntimeProperties::GENE_ID_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Ident: {} Gene ident file not defined", nomenclature_ident);

  }


  if (nomenclature_ident == RuntimeProperties::NOMENCLATURE_UNIPROTID) {

    ParseUniprotId parse_uniprot;
    if (not parse_uniprot.parseUniprotFile(file_name_opt.value())) {

      ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Unable to parse Uniprot Nomenclature file: {}", file_name_opt.value());

    }

    std::shared_ptr<const UniprotResource> gene_id_resource(std::make_shared<const UniprotResource>(nomenclature_ident, parse_uniprot.getUniproResource()));

    resource_ptr->addResource(gene_id_resource);

  } else if (nomenclature_ident == RuntimeProperties::NOMENCLATURE_ENSEMBL) {

    ParseGeneIdents parse_gene_idents;
    if (not parse_gene_idents.parseIdentFile(file_name_opt.value())) {

      ExecEnv::log().critical("ExecutePackage::loadGeneNomenclatureResource, Unable to parse Ensembl Nomenclature file: {}", file_name_opt.value());

    }

    std::shared_ptr<const EnsemblHGNCResource> gene_id_resource(std::make_shared<const EnsemblHGNCResource>(nomenclature_ident, parse_gene_idents.getSynonymVector()));

    resource_ptr->addResource(gene_id_resource);

  } else {

    ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Nomenclature resource ident: {} unknown", nomenclature_ident);

  }

}


void kgl::ExecutePackage::loadHsGenomeGenealogyResource(const std::string& resource_type,
                                                        const std::string& genealogy_ident,
                                                        const std::shared_ptr<AnalysisResources>& resource_ptr) const {


  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, genealogy_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadHsGenomeGenealogyResource, Genome Genealogy Database: {}, not defined", genealogy_ident);

  }

  auto const& params = params_opt.value();
  auto file_name_opt = params.getParameter(RuntimeProperties::GENEALOGY_ID_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Ident: {} Gene ident file not defined", genealogy_ident);

  }

  std::shared_ptr<HsGenomeGenealogyData> genealogy_data(std::make_shared<HsGenomeGenealogyData>(genealogy_ident));
  ParseHsGenomeGenealogyFile genealogy_parser(genealogy_data);
  genealogy_parser.readParseImpl(file_name_opt.value());

  resource_ptr->addResource(genealogy_data);

}


void kgl::ExecutePackage::loadHsGenomeAuxResource(const std::string& resource_type, const std::string& genome_aux_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

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


void kgl::ExecutePackage::loadCitationResource(const std::string& resource_type, const std::string& citation_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

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


void kgl::ExecutePackage::loadEntrezGeneResource(const std::string& resource_type, const std::string& entrez_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

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



void kgl::ExecutePackage::loadPubmedAPIResource(const std::string& resource_type, const std::string& api_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(api_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadPubmedAPIResource, Pubmed API Database: {}, not defined", api_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto api_resource_ptr = std::dynamic_pointer_cast<const RuntimePubmedAPIResource>(resource_base_ptr);

  if (not api_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadPubmedAPIResource, Resource: {} is not a Pubmed API Resource", resource_ident);

  }

  std::shared_ptr<PubmedRequester> pubmed_api_ptr(std::make_shared<PubmedRequester>( api_resource_ptr->apiIdentifier(),
                                                                                     api_resource_ptr->publicationCacheFile(),
                                                                                     api_resource_ptr->citationCacheFile()));

  resource_ptr->addResource(pubmed_api_ptr);

}


void kgl::ExecutePackage::loadPf7SampleResource(const std::string& resource_type, const std::string& Pf7_sample_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto result = runtime_config_.resourceMap().find(Pf7_sample_ident);
  if (result == runtime_config_.resourceMap().end()) {

    ExecEnv::log().critical("ExecutePackage::loadPf7SampleResource, Pf7 Sample Data: {}, not defined", Pf7_sample_ident);

  }

  auto const& [resource_ident, resource_base_ptr] = *result;
  auto Pf7_resource_ptr = std::dynamic_pointer_cast<const RuntimePf7Resource>(resource_base_ptr);

  if (not Pf7_resource_ptr) {

    ExecEnv::log().critical("ExecutePackage::loadPf7SampleResource, Resource: {} is not a  f7 Sample Data Resource", resource_ident);

  }

  ParsePf7Sample Pf7_sample_parser;
  if (not Pf7_sample_parser.parsePf7SampleFile(Pf7_resource_ptr->Pf7FileName())) {

    ExecEnv::log().critical("ExecutePackage::loadPf7SampleResource; failed to create Entrez resource from file: {}", Pf7_resource_ptr->Pf7FileName());

  }

  std::shared_ptr<Pf7SampleResource> Pf7Sample_ptr(std::make_shared<Pf7SampleResource>(Pf7_resource_ptr->Pf7Identifier(), Pf7_sample_parser.getPf7SampleVector()));

  resource_ptr->addResource(Pf7Sample_ptr);

}

