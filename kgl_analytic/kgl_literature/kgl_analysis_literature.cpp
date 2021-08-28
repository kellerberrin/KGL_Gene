//
// Created by kellerberrin on 24/8/21.
//

#include "kgl_mutation/kgl_analysis_mutation_data.h"
#include "kgl_analysis_literature.h"
#include "kgl_analysis_gene_sequence.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::LiteratureAnalysis::initializeAnalysis(const std::string& work_directory,
                                                 const ActiveParameterList& named_parameters,
                                                 const std::shared_ptr<const AnalysisResources>& resource_ptr) {


  // Setup and clear the directories to hold analysis output.
  // Most output in here.
  ident_work_directory_ = work_directory + std::string("/") + ident();
  if (not Utility::createDirectory(ident_work_directory_)) {

    ExecEnv::log().critical("MutationAnalysis::initializeAnalysis, unable to create analysis results directory: {}",
                            ident_work_directory_);

  }
  // Gene specific output in here.
  ident_gene_work_directory_ = ident_work_directory_ + std::string("/") + gene_output_directory;
  if (not Utility::recreateDirectory(ident_gene_work_directory_)) {

    ExecEnv::log().critical("MutationAnalysis::initializeAnalysis, unable to recreate gene analysis results directory: {}",
                            ident_gene_work_directory_);

  }


  // Get the analysis parameters.
  ExecEnv::log().info("Default Analysis Id: {} initialized with analysis results directory: {}", ident(), ident_work_directory_);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  // Get the analysis resources.
  ref_genome_ptr_ = resource_ptr->getSingleResource<const GenomeReference>(RuntimeResourceType::GENOME_DATABASE);
  uniprot_nomenclature_ptr_ = resource_ptr->getSingleResource<const UniprotResource>(RuntimeResourceType::GENE_NOMENCLATURE, ResourceBase::NOMENCLATURE_UNIPROTID);
  ensembl_nomenclature_ptr_ = resource_ptr->getSingleResource<const EnsemblHGNCResource>(RuntimeResourceType::GENE_NOMENCLATURE, ResourceBase::NOMENCLATURE_ENSEMBL);
  entrez_nomenclature_ptr_ = resource_ptr->getSingleResource<const EntrezResource>(RuntimeResourceType::ENTREZ_GENE);
  allele_citation_ptr_ = resource_ptr->getSingleResource<const CitationResource>(RuntimeResourceType::ALLELE_CITATION);
  pubmed_requestor_ptr_ = resource_ptr->getSingleResource<const PubmedRequester>(RuntimeResourceType::PUBMED_API);

  // Create the gene list.
  gene_literature_.defineGenes( ref_genome_ptr_, uniprot_nomenclature_ptr_, entrez_nomenclature_ptr_);

  return true;

}



// Perform pre-processing (generally just type casting) for each file read into the analysis object.
bool kgl::LiteratureAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

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
bool kgl::LiteratureAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());


  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::LiteratureAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  const size_t minimum_publication_count{10};
  gene_literature_.outputGenePmid(pubmed_requestor_ptr_, ident_gene_work_directory_, minimum_publication_count);

  const size_t max_genes{1000000};
  const size_t min_genes{0};
  const size_t min_citations{0};
  gene_literature_.outputPmidGene(pubmed_requestor_ptr_, ident_work_directory_, max_genes, min_genes, min_citations);

  return true;

}

