//
// Created by kellerberrin on 24/8/21.
//

#include "kga_analysis_mutation_data.h"
#include "kga_analysis_literature.h"
#include "kgl_analysis_gene_sequence.h"
#include "kga_analysis_literature_publication.h"


namespace kga = kellerberrin::genome::analysis;


// Setup the analytics to process VCF data.
bool kga::LiteratureAnalysis::initializeAnalysis(const std::string& work_directory,
                                                 const ActiveParameterList& named_parameters,
                                                 const std::shared_ptr<const AnalysisResources>& resource_ptr) {


  // Setup and clear the directories to hold analysis output.
  // The top level directory for this analysis type.
  // This directory and sub-directories are recreated each time the analysis is executed.
  // ALL PREVIOUS ANALYSIS FILES ARE DELETED.
  ident_work_directory_ = work_directory + std::string("/") + ident();
  if (not Utility::createDirectory(ident_work_directory_)) {

    ExecEnv::log().critical("MutationAnalysis::initializeAnalysis, unable to create analysis results directory: {}",
                            ident_work_directory_);

  }

  // Gene specific output in here.
  gene_work_directory_ = ident_work_directory_ + std::string("/") + GENE_DIRECTORY_NAME_;
  if (not Utility::recreateDirectory(gene_work_directory_)) {

    ExecEnv::log().critical("MutationAnalysis::initializeAnalysis, unable to recreate gene analysis results directory: {}",
                            gene_work_directory_);

  }

  // Literature citation, journal, author statistics output in here.
  analysis_work_directory_ = ident_work_directory_ + std::string("/") + ANALYSIS_DIRECTORY_NAME_;
  if (not Utility::createDirectory(analysis_work_directory_)) {

    ExecEnv::log().critical("MutationAnalysis::initializeAnalysis, unable to recreate gene analysis results directory: {}",
                            analysis_work_directory_);

  }

  // Get the analysis parameters.
  ExecEnv::log().info("Default Analysis Id: {} initialized with analysis results directory: {}", ident(), ident_work_directory_);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  // Get the analysis resources.
  ref_genome_ptr_ = resource_ptr->getSingleResource<const GenomeReference>(ResourceProperties::GENOME_RESOURCE_ID_);
  uniprot_nomenclature_ptr_ = resource_ptr->getSingleResource<const UniprotResource>(ResourceProperties::GENE_NOMENCLATURE_RESOURCE_ID_,
                                                                                     ResourceProperties::NOMENCLATURE_UNIPROTID);
  ensembl_nomenclature_ptr_ = resource_ptr->getSingleResource<const EnsemblHGNCResource>(ResourceProperties::GENE_NOMENCLATURE_RESOURCE_ID_,
                                                                                         ResourceProperties::NOMENCLATURE_ENSEMBL);
  entrez_nomenclature_ptr_ = resource_ptr->getSingleResource<const EntrezResource>(ResourceProperties::ENTREZ_RESOURCE_ID_);
  allele_citation_ptr_ = resource_ptr->getSingleResource<const CitationResource>(ResourceProperties::CITATION_RESOURCE_ID_);
  pubmed_requestor_ptr_ = resource_ptr->getSingleResource<const PubmedRequester>(ResourceProperties::PUBMED_API_RESOURCE_ID_);

  // Create the gene list.
  gene_literature_.defineGenes( ref_genome_ptr_, uniprot_nomenclature_ptr_, entrez_nomenclature_ptr_);

  return true;

}



// Perform pre-processing (generally just type casting) for each file read into the analysis object.
bool kga::LiteratureAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called with file: {}", ident(), data_ptr->fileId());

  auto file_characteristic = data_ptr->dataCharacteristic();

  if (file_characteristic.data_source == DataSourceEnum::BioPMID) {

    auto bio_pmid_ptr = std::dynamic_pointer_cast<const BioPMIDFileData>(data_ptr);

    if (not bio_pmid_ptr) {

      ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast data file to Bio PMID object, severe error.");

    }

    ExecEnv::log().info("Bio PMID MeSH Disease Map Size: {}", bio_pmid_ptr->diseaseMeSHMap().size());
    ExecEnv::log().info("Bio PMID Entrez Gene Map Size: {}", bio_pmid_ptr->entrezMap().size());

    auto const disease_pmid_set = bio_pmid_ptr->selectDiseaseBioPMID(MutationAnalysisData::malariaMeSHList());
    gene_literature_.updatePMIDStatistics(disease_pmid_set, bio_pmid_ptr);

  } else {

    ExecEnv::log().error("MutationAnalysis::fileReadAnalysis unknown file type: {}", file_characteristic.source_text);

  }

  return true;

}


// Perform the genetic analysis per iteration.
bool kga::LiteratureAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All file data has been presented, finalize analysis and write results.
bool kga::LiteratureAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  const size_t minimum_publication_count{10};
  // Output the publications relevant to each gene.
  gene_literature_.outputGenePmid(pubmed_requestor_ptr_, gene_work_directory_, minimum_publication_count);

  const size_t max_genes{1000000};
  const size_t min_genes{0};
  const size_t min_citations{0};
  // Output the Genes relevant for each publication.
  gene_literature_.outputPmidGene(pubmed_requestor_ptr_, ident_work_directory_, max_genes, min_genes, min_citations);

  // Output literature analytics
  PublicationLiterature publication_literature(pubmed_requestor_ptr_->getAllCachedPublications());
  publication_literature.writeAuthorAnalysis(analysis_work_directory_);
  publication_literature.writeYearAnalysis(analysis_work_directory_);
  publication_literature.writeJournalAnalysis(analysis_work_directory_);
  publication_literature.writeCitationPeriod(analysis_work_directory_);
  publication_literature.writeCitationVariance(analysis_work_directory_);
  publication_literature.writeCitationQuantiles(analysis_work_directory_);
  publication_literature.writeCitationHistogram(analysis_work_directory_);
  publication_literature.writeCitationData(analysis_work_directory_);

  publication_literature.writePublicationCitations(analysis_work_directory_, "22157630");
  publication_literature.writePublicationCitations(analysis_work_directory_, "18717995");
  publication_literature.writePublicationCitations(analysis_work_directory_, "7886766");

  auto publication_ptr = publication_literature.mostRecentPublication();
  ExecEnv::log().info("LiteratureAnalysis::finalizeAnalysis; most recent publication pmid: {} was on: {}",
                      publication_ptr->pmid(), publication_ptr->publicationDate().text());

  return true;

}

