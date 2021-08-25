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


  // Get the analysis parameters.
  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
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
  // For efficiency, the Pubmed API requestor uses file caches and needs to know the runtime directory structure.
  pubmed_requestor_ptr_->setWorkDirectory(work_directory);
  work_directory_ = work_directory;
  literature_directory_ = work_directory + LITERATURE_DIRECTORY_;
  // Create the literature directory if it does not exist.
  if (not Utility::createDirectory(literature_directory_)) {

    ExecEnv::log().critical( "LiteratureAnalysis::initializeAnalysis, unable to create literature directory: {}",
                             literature_directory_);

  }

  defineGenes( ref_genome_ptr_, uniprot_nomenclature_ptr_, entrez_nomenclature_ptr_);

  return true;

}




// Perform the genetic analysis per iteration.
void kgl::LiteratureAnalysis::defineGenes( const std::shared_ptr<const GenomeReference>& genome_ptr,
                                           const std::shared_ptr<const UniprotResource>& uniprot_nomenclature_ptr,
                                           const std::shared_ptr<const EntrezResource>& entrez_nomenclature_ptr) {

  // Only execute this function once.


  gene_vector_.clear();

  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [offset, gene_ptr] : contig_ptr->getGeneMap()) {

      std::vector<std::string> name_vec;
      gene_ptr->getAttributes().getName(name_vec);
      std::string symbol_id;
      std::string gaf_id;

      if (not name_vec.empty()) {

        symbol_id = name_vec.front();

      }

      std::string hgnc_id = gene_ptr->getAttributes().getHGNC();
      std::vector<std::string> ensembl_ids = uniprot_nomenclature_ptr->HGNCToEnsembl(hgnc_id);

      std::string entrez_id = entrez_nomenclature_ptr->symbolToEntrez(symbol_id);

      GeneCharacteristic gene_record;
      gene_record.geneDefinition(gene_ptr, genome_ptr->genomeId(), symbol_id, hgnc_id, ensembl_ids, gaf_id, entrez_id);
      gene_vector_.push_back(gene_record);

    } // Gene.

  } // Contig.

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
    updatePMIDStatistics(disease_pmid_set, bio_pmid_ptr);

  } else {

    ExecEnv::log().error("MutationAnalysis::fileReadAnalysis unknown file type: {}", file_characteristic.source_text);

  }

  return true;

}


void kgl::LiteratureAnalysis::updatePMIDStatistics( const std::set<std::string>& disease_pmid_set,
                                                    const std::shared_ptr<const BioPMIDFileData>& bio_pmid_ptr) {

  size_t gene_disease_count{0};
  size_t entrez_empty{0};
  for (auto& gene_data : gene_vector_) {

    const std::string& entrez_id = gene_data.entrezId();

    if (entrez_id.empty()) {

      ++entrez_empty;
      continue;

    }

    auto const entrez_pmid = bio_pmid_ptr->entrezPMID(entrez_id);
    std::set<std::string> gene_disease_pmids;
    for (auto const& pmid : entrez_pmid) {

      if (disease_pmid_set.contains(pmid)) {

        gene_disease_pmids.insert(pmid);

      }

    }

    gene_disease_count += gene_disease_pmids.size();
    gene_data.update_pmid(entrez_pmid.size(), gene_disease_pmids);

  }

  ExecEnv::log().info("LiteratureAnalysis::updatePMIDStatistics; Gene with Empty Entrez Id: {}, Gene Disease pmid: {}", entrez_empty, gene_disease_count);

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
  outputGenePmid(literature_directory_, minimum_publication_count);

  const size_t max_genes{1000000};
  const size_t min_genes{0};
  const size_t min_citations{10};
  outputPmidGene(literature_directory_, max_genes, min_genes, min_citations);

  return true;

}

// Ranks the publications by number of gene references.
void kgl::LiteratureAnalysis::outputGenePmid(const std::string& literature_directory, size_t pmid_count) const {

  size_t gene_lit_count{0};
  size_t gene_lit_total{0};
  for (auto const& gene : gene_vector_) {

    if (not gene.diseasePublications().empty()) {

      ++gene_lit_total;

      if (gene.diseasePublications().size() >= pmid_count) {

        std::string literature_file = std::string(gene.symbolId()) + std::string(".txt");
        literature_file = Utility::filePath(literature_file, literature_directory);
        std::ofstream out_file(literature_file);


        if (out_file.good()) {

          ++gene_lit_count;
          gene.writeGenePublications(out_file, pubmed_requestor_ptr_);

        } else {

          ExecEnv::log().error("LiteratureAnalysis::outputGenePmid; problem opening file: {}", literature_file);

        }

      } // if publications >= pmid_count

    } // if publication

  } // for all genes.

  ExecEnv::log().info("Literature Analysis; Total Genes: {}, Total Lit Gene: {}, Lit Gene Files: {}",
                      gene_vector_.size(), gene_lit_total, gene_lit_count);

}

// Ranks the publications by number of gene references.
void kgl::LiteratureAnalysis::outputPmidGene(const std::string& literature_directory,
                                             size_t max_genes,
                                             size_t min_genes,
                                             size_t min_citations) const {

  std::map<std::string, std::set<std::string>> reference_map;
  for (auto const& gene : gene_vector_) {

    for (auto const& pmid : gene.diseasePublications()) {

      auto result = reference_map.find(pmid);
      if (result == reference_map.end()) {

        reference_map.emplace(pmid, std::set<std::string>{gene.symbolId()});

      } else {

        auto& [pmid_key, gene_set] = *result;
        gene_set.insert(gene.symbolId());

      }

    }

  }


  std::vector<std::string> pmid_vector;
  for (auto const& [pmid, gene_set] : reference_map) {

    pmid_vector.push_back(pmid);

  }

  // Get all the publications.
  auto publication_map = pubmed_requestor_ptr_->getCachedPublications(pmid_vector);

  std::multimap<size_t, std::string> count_map;
  size_t fail_filter{0};
  size_t pass_filter{0};
  for (auto const& [pmid, gene_set] : reference_map) {

    auto pub_result = publication_map.find(pmid);
    if (pub_result == publication_map.end()) {

      ExecEnv::log().error("LiteratureAnalysis::outputPmidGene; unable to find pmid publication: {}", pmid);
      continue;

    }

    auto const& [pub_pmid, publication] = *pub_result;

    if (publication.citedBy().size() >= min_citations) {

      if (gene_set.size() <= max_genes and gene_set.size() >= min_genes) {


        bool filter_result = filterPublication(publication);

        if (filter_result) {

          ++pass_filter;
          count_map.emplace(gene_set.size(), pmid);

        } else {

          ++fail_filter;

        }

      } // if gene count

    } // if citations

  }

  ExecEnv::log().info( "LiteratureAnalysis::outputPmidGene; pass filter: {}, fail filter: {}", pass_filter, fail_filter);

  std::string literature_file{"reference_count.txt"};
  literature_file = Utility::filePath(literature_file, literature_directory);
  std::ofstream out_file(literature_file);

  if (out_file.good()) {

    // print the most referenced publications first.
    for (auto iter = count_map.rbegin(); iter != count_map.rend(); ++iter) {

      auto const& [count, pmid] = *iter;

      auto ref_result = reference_map.find(pmid);
      if (ref_result == reference_map.end()) {

        ExecEnv::log().error("LiteratureAnalysis::outputPmidGene; unable to find pmid reference: {}", pmid);
        continue;

      }

      auto pub_result = publication_map.find(pmid);
      if (pub_result == publication_map.end()) {

        ExecEnv::log().error("LiteratureAnalysis::outputPmidGene; unable to find pmid publication: {}", pmid);
        continue;

      }

      auto const& [pub_pmid, publication] = *pub_result;
      auto const& [ref_pmid, gene_set] = *ref_result;

      out_file << "*************************************************************\n";
      out_file << "Gene count: " << count << '\n';
      const size_t genes_per_line{20};
      size_t gene_count{0};
      for (auto const& gene : gene_set) {

        out_file << gene;
        if (gene != *gene_set.rbegin()) {

          out_file << ", ";

        }
        ++gene_count;
        if (gene_count % genes_per_line == 0) {

          out_file << '\n';

        }

      }

      out_file << "\n\n";

      publication.output(out_file);

      out_file << "\n*************************************************************\n";

    }

  } else {

    ExecEnv::log().error("LiteratureAnalysis::outputPmidGene; problem opening file: {}", literature_file);

  }

}

bool kgl::LiteratureAnalysis::filterPublication(const PubMedPublicationSummary& publication) const {

  static const std::vector<std::string> Mesh_codes{ "D010963" /* Plasmodium falciparum */
                                                   ,"D016778" /*  Malaria, Falciparum */
                                                   , "D008288" /* Malaria */ };

  static const std::vector<std::string> search_text { "Plasmodium", "plasmodium", "Falciparum",  "falciparum", "Malaria", "malaria" };

  if (publication.hasMeSHCode(Mesh_codes)) {

    return true;

  }

  if (publication.hasTitleText(search_text)) {

    return true;

  }

  if (publication.hasAbstractText(search_text)) {

    return true;

  }

  return false;

}
