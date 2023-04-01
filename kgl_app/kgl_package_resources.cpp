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

void kgl::ExecutePackage::loadHsGenomeAuxResource(const std::string& resource_type,
                                                  const std::string& genome_aux_ident,
                                                  const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, genome_aux_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadHsGenomeAuxResource, Genome Genealogy Database: {}, not defined", genome_aux_ident);

  }

  auto const& params = params_opt.value();
  auto file_name_opt = params.getParameter(RuntimeProperties::AUX_ID_FILE_);
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
  auto file_name_opt = params.getParameter(RuntimeProperties::CITATION_FILE_);
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
  auto file_name_opt = params.getParameter(RuntimeProperties::ENTREZ_FILE_);
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

  auto publication_file_opt = params.getParameter(RuntimeProperties::PUBMED_PUBLICATION_CACHE_);
  if (not publication_file_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPubmedAPIResource, Ident: {} Publication file not defined", api_ident);

  }

  auto citation_file_opt = params.getParameter(RuntimeProperties::PUBMED_CITATION_CACHE_);
  if (not citation_file_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPubmedAPIResource, Ident: {} Citation file not defined", api_ident);

  }

  std::shared_ptr<PubmedRequester> pubmed_api_ptr(std::make_shared<PubmedRequester>( api_ident,
                                                                                     publication_file_opt.value(),
                                                                                     citation_file_opt.value()));

  resource_ptr->addResource(pubmed_api_ptr);

}

void kgl::ExecutePackage::loadPf7SampleResource(const std::string& resource_type,
                                                const std::string& Pf7_sample_ident,
                                                const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, Pf7_sample_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPf7SampleResource, Pf7 Sample Data: {}, not defined", Pf7_sample_ident);

  }

  auto const& params = params_opt.value();
  auto file_name_opt = params.getParameter(RuntimeProperties::PF7_SAMPLE_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPf7SampleResource, Ident: {} Pf7 Sample Data file not defined", Pf7_sample_ident);

  }

  ParsePf7Sample Pf7_sample_parser;
  if (not Pf7_sample_parser.parsePf7SampleFile(file_name_opt.value())) {

    ExecEnv::log().critical("ExecutePackage::loadPf7SampleResource; failed to create Pf7 Sample Data resource from file: {}", file_name_opt.value());

  }

  std::shared_ptr<Pf7SampleResource> Pf7Sample_ptr(std::make_shared<Pf7SampleResource>(Pf7_sample_ident, Pf7_sample_parser.getPf7SampleVector()));

  resource_ptr->addResource(Pf7Sample_ptr);

}

