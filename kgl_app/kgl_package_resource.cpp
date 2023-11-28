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
#include "kgl_hsgenealogy_parser.h"
#include "kgl_pubmed_resource.h"


namespace kgl = kellerberrin::genome;


std::shared_ptr<const kgl::AnalysisResources> kgl::ExecutePackage::loadRuntimeResources(const RuntimePackage& package) const {

  std::shared_ptr<AnalysisResources> resource_ptr(std::make_shared<AnalysisResources>());

  for (auto const& [resource_type, resource_ident]  : package.resourceList()) {

    switch (Utility::hash(resource_type)) {

      case Utility::hash(ResourceProperties::GENOME_RESOURCE_ID_):
        loadGenomeResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::ONTOLOGY_RESOURCE_ID_):
        loadOntologyResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::GENE_NOMENCLATURE_RESOURCE_ID_):
        loadHsGeneNomenclatureResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::GENEALOGY_RESOURCE_ID_):
        loadHsGenomeGenealogyResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::GENOMEAUX_RESOURCE_ID_):
        loadHsGenomeAuxResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::CITATION_RESOURCE_ID_):
        loadCitationResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::ENTREZ_RESOURCE_ID_):
        loadEntrezGeneResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::PUBMED_API_RESOURCE_ID_):
        loadPubmedAPIResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::PF7SAMPLE_RESOURCE_ID_):
        loadPf7SampleResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::PF7FWS_RESOURCE_ID_):
        loadPf7FwsResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::PF7DISTANCE_RESOURCE_ID_):
        loadPf7DistanceResource(resource_type, resource_ident, resource_ptr);
        break;

      case Utility::hash(ResourceProperties::PF3K_COI_RESOURCE_ID_):
        loadPf3KCOIResource(resource_type, resource_ident, resource_ptr);
        break;

      default:
        ExecEnv::log().critical("ExecutePackage::loadRuntimeResources, Package: {} Attempt to load unknown resource type: '{}', resource ident: '{}'.",
                                package.packageIdentifier(), resource_type, resource_ident);
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
  auto fasta_opt = params.getParameter(ResourceProperties::FASTA_FILE_);
  if (not fasta_opt) {

    ExecEnv::log().critical("ExecutePackage::loadGenomeResource, Reference Genome: {}, Fasta file not defined", genome_ident);

  }
  auto gff_opt = params.getParameter(ResourceProperties::GFF_FILE_);
  if (not gff_opt) {

    ExecEnv::log().critical("ExecutePackage::loadGenomeResource, Reference Genome: {}, Gff file not defined", genome_ident);

  }
  auto translation_opt = params.getParameter(ResourceProperties::TRANSLATION_TABLE_);
  if (not translation_opt) {

    ExecEnv::log().critical("ExecutePackage::loadGenomeResource, Reference Genome: {}, Translation table not defined", genome_ident);

  }
  // The Gaf file parameter is optional
  auto gaf_opt = params.getParameter(ResourceProperties::GAF_ANNOTATION_FILE_);
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
  auto go_graph_opt = params.getParameter(ResourceProperties::GO_ONTOLOGY_FILE_);
  if (not go_graph_opt) {

    ExecEnv::log().critical("ExecutePackage::loadOntologyResource, Ontology Ident: {} Go graph file not defined", ontology_ident);

  }

  auto gaf_opt = params.getParameter(ResourceProperties::GAF_ANNOTATION_FILE_);
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
  auto file_name_opt = params.getParameter(ResourceProperties::GENE_NOMENCLATURE_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Ident: {} Gene ident file not defined", nomenclature_ident);

  }

  if (nomenclature_ident == ResourceProperties::NOMENCLATURE_UNIPROTID) {

    ParseUniprotId parse_uniprot;
    if (not parse_uniprot.parseUniprotFile(file_name_opt.value())) {

      ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Unable to parse Uniprot Nomenclature file: {}", file_name_opt.value());

    }

    std::shared_ptr<const UniprotResource> gene_id_resource(std::make_shared<const UniprotResource>(nomenclature_ident, parse_uniprot.getUniproResource()));

    resource_ptr->addResource(gene_id_resource);

  } else if (nomenclature_ident == ResourceProperties::NOMENCLATURE_ENSEMBL) {

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
  auto file_name_opt = params.getParameter(ResourceProperties::GENEALOGY_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("ExecutePackage::loadHsGeneNomenclatureResource, Ident: {} Gene ident file not defined", genealogy_ident);

  }

  std::shared_ptr<HsGenomeGenealogyData> genealogy_data(std::make_shared<HsGenomeGenealogyData>(genealogy_ident));
  ParseHsGenomeGenealogyFile genealogy_parser(genealogy_data);
  genealogy_parser.readParseImpl(file_name_opt.value());

  resource_ptr->addResource(genealogy_data);

}

void kgl::ExecutePackage::loadHsGenomeAuxResource(const std::string& resource_type,
                                                  const std::string& genome_aux_ident,
                                                  const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, genome_aux_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadHsGenomeAuxResource, Genome Genealogy Database: {}, not defined", genome_aux_ident);

  }

  auto const& params = params_opt.value();
  auto file_name_opt = params.getParameter(ResourceProperties::GENOMEAUX_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("ExecutePackage::loadHsGenomeAuxResource, Ident: {} Aux file not defined", genome_aux_ident);

  }

  std::shared_ptr<HsGenomeAuxData> genome_aux_data(std::make_shared<HsGenomeAuxData>(genome_aux_ident));
  ParseHsGenomeAuxFile genome_aux_parser(genome_aux_data);
  genome_aux_parser.readParseImpl(file_name_opt.value());

  resource_ptr->addResource(genome_aux_data);

}

void kgl::ExecutePackage::loadCitationResource(const std::string& resource_type,
                                               const std::string& citation_ident,
                                               const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, citation_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadCitationResource, Allele Citation Database: {}, not defined", citation_ident);

  }

  auto const& params = params_opt.value();
  auto file_name_opt = params.getParameter(ResourceProperties::CITATION_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("ExecutePackage::loadCitationResource, Ident: {} Citation file not defined", citation_ident);

  }

  ParseCitations citation_parser;
  if (not citation_parser.parseCitationFile(file_name_opt.value())) {

    ExecEnv::log().critical("ExecutePackage::loadCitationResource; failed to create citation resource from file: {}", file_name_opt.value());

  }

  std::shared_ptr<CitationResource> citation_ptr(std::make_shared<CitationResource>(citation_ident, citation_parser.getCitationMap()));

  resource_ptr->addResource(citation_ptr);

}

void kgl::ExecutePackage::loadEntrezGeneResource(const std::string& resource_type,
                                                 const std::string& entrez_ident,
                                                 const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, entrez_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadEntrezGeneResource, Entrez Gene Database: {}, not defined", entrez_ident);

  }

  auto const& params = params_opt.value();
  auto file_name_opt = params.getParameter(ResourceProperties::ENTREZ_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("EExecutePackage::loadEntrezGeneResource, Ident: {} Entrez file not defined", entrez_ident);

  }

  ParseEntrez entrez_parser;
  if (not entrez_parser.parseEntrezFile(file_name_opt.value())) {

    ExecEnv::log().critical("ExecutePackage::loadEntrezGeneResource; failed to create Entrez resource from file: {}", file_name_opt.value());

  }

  std::shared_ptr<EntrezResource> entrez_ptr(std::make_shared<EntrezResource>(entrez_ident, entrez_parser.getEntrezVector()));

  resource_ptr->addResource(entrez_ptr);

}

void kgl::ExecutePackage::loadPubmedAPIResource(const std::string& resource_type,
                                                const std::string& api_ident,
                                                const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, api_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPubmedAPIResource, Pubmed API Database: {}, not defined", api_ident);

  }

  auto const& params = params_opt.value();

  auto publication_file_opt = params.getParameter(ResourceProperties::PUBMED_PUBLICATION_CACHE_);
  if (not publication_file_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPubmedAPIResource, Ident: {} Publication file not defined", api_ident);

  }

  auto citation_file_opt = params.getParameter(ResourceProperties::PUBMED_CITATION_CACHE_);
  if (not citation_file_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPubmedAPIResource, Ident: {} Citation file not defined", api_ident);

  }

  std::shared_ptr<PubmedRequester> pubmed_api_ptr(std::make_shared<PubmedRequester>( api_ident,
                                                                                     publication_file_opt.value(),
                                                                                     citation_file_opt.value()));

  resource_ptr->addResource(pubmed_api_ptr);

}
