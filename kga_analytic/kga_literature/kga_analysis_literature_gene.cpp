//
// Created by kellerberrin on 27/8/21.
//

#include "kga_analysis_literature_gene.h"


namespace kga = kellerberrin::genome::analysis;



// Perform the genetic analysis per iteration.
void kga::GeneLiterature::defineGenes( const std::shared_ptr<const GenomeReference>& genome_ptr,
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



void kga::GeneLiterature::updatePMIDStatistics( const std::set<std::string>& disease_pmid_set,
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
    gene_data.update_pmid(entrez_pmid, gene_disease_pmids);

  }

  ExecEnv::log().info("LiteratureAnalysis::updatePMIDStatistics; Gene with Empty Entrez Id: {}, Gene Disease pmid: {}", entrez_empty, gene_disease_count);

}


// Ranks the publications by number of gene references.
void kga::GeneLiterature::outputGenePmid(const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr,
                                         const std::string& literature_directory,
                                         size_t pmid_count) const {

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
          writeGenePublications(out_file, gene, pubmed_requestor_ptr);

        } else {

          ExecEnv::log().error("LiteratureAnalysis::outputGenePmid; problem opening file: {}", literature_file);

        }

      } // if publications >= pmid_count

    } // if publication

  } // for all genes.

  ExecEnv::log().info("Literature Analysis; Total Genes: {}, Total Lit Gene: {}, Lit Gene Files: {}",
                      gene_vector_.size(), gene_lit_total, gene_lit_count);

}


void kga::GeneLiterature::writeGenePublications( std::ostream& out_file,
                                                 const GeneCharacteristic& gene,
                                                 const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr) {


  const char output_delimiter{'\n'};
  const char CONCAT_TOKEN_{'&'};

  out_file << gene.genome() << output_delimiter
  << gene.contigId() << output_delimiter
  << gene.entrezId() << output_delimiter
  << gene.HGNC() << output_delimiter;

  std::string ensembl_text;
  for (auto const& ensembl_id : gene.ensemblIds()) {

    ensembl_text += ensembl_id;
    if (ensembl_id != gene.ensemblIds().back()) {

      ensembl_text += CONCAT_TOKEN_;

    }

  }

  out_file << ensembl_text << output_delimiter
  << gene.gafId() << output_delimiter
  << gene.symbolId() << output_delimiter
  << gene.description() << output_delimiter
  << gene.citations() << output_delimiter
  << gene.diseaseCites().size() << output_delimiter;

  out_file << output_delimiter
  << "************************************************************************************";

  std::vector<std::string> pmid_vector;
  for (auto const& pmid  :  gene.diseaseCites()) {

    pmid_vector.push_back(pmid);

  }

  auto literature_map = pubmed_requestor_ptr->getCachedPublications(pmid_vector);

  // Resort the literature map by number of citations.
  std::multimap<size_t, std::shared_ptr<const PublicationSummary>> citation_rank_map;
  for (auto const& [pmid, publication_ptr] : literature_map) {

    if (filterPublication(*publication_ptr)) {

      citation_rank_map.emplace(publication_ptr->citedBy().size(), publication_ptr);

    }

  }

  // Print most cited articles first.
  for (auto iter = citation_rank_map.rbegin(); iter != citation_rank_map.rend(); ++iter) {

    auto const& [cites, publication_ptr] = *iter;

    out_file << "\n******************************************" << output_delimiter;

    publication_ptr->extendedBiblio(out_file);

    out_file << "\n******************************************" << output_delimiter;

  }

}


// Ranks the publications by number of gene references.
void kga::GeneLiterature::outputPmidGene( const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr,
                                          const std::string& literature_directory,
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
  auto publication_map = pubmed_requestor_ptr->getCachedPublications(pmid_vector);

  std::multimap<size_t, std::string> count_map;
  size_t fail_filter{0};
  size_t pass_filter{0};
  for (auto const& [pmid, gene_set] : reference_map) {

    auto pub_result = publication_map.find(pmid);
    if (pub_result == publication_map.end()) {

      ExecEnv::log().error("LiteratureAnalysis::outputPmidGene; unable to find pmid publication_ptr: {}", pmid);
      continue;

    }

    auto const& [pub_pmid, publication_ptr] = *pub_result;

    if (publication_ptr->citedBy().size() >= min_citations) {

      if (gene_set.size() <= max_genes and gene_set.size() >= min_genes) {


        bool filter_result = filterPublication(*publication_ptr);

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

        ExecEnv::log().error("LiteratureAnalysis::outputPmidGene; unable to find pmid publication_ptr: {}", pmid);
        continue;

      }

      auto const& [pub_pmid, publication_ptr] = *pub_result;
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

      publication_ptr->extendedBiblio(out_file);

      out_file << "\n*************************************************************\n";

    }

  } else {

    ExecEnv::log().error("LiteratureAnalysis::outputPmidGene; problem opening file: {}", literature_file);

  }

}

