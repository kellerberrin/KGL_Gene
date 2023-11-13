//
// Created by kellerberrin on 31/8/21.
//

#include "kga_analysis_mutation_gene_allele_pop.h"
#include "kgl_literature_filter.h"


namespace kga = kellerberrin::genome::analysis;


void kga::GeneratePopulationAllele::initialize( const std::shared_ptr<const HsGenomeAux>& genome_aux_ptr,
                                                const std::shared_ptr<const UniprotResource>& uniprot_nomenclature_ptr,
                                                const std::shared_ptr<const EntrezResource>& entrez_nomenclature_ptr,
                                                const std::shared_ptr<const CitationResource>& allele_citation_ptr,
                                                const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr) {

  // Initialize all the resource pointers.
  genome_aux_ptr_ = genome_aux_ptr;
  uniprot_nomenclature_ptr_ = uniprot_nomenclature_ptr;
  entrez_nomenclature_ptr_ = entrez_nomenclature_ptr;
  allele_citation_ptr_ = allele_citation_ptr;
  pubmed_requestor_ptr_ = pubmed_requestor_ptr;

  // Initialize the reference ethnic analysis object.
  reference_ethnic_.updatePopulations(genome_aux_ptr_);

}


void kga::GeneratePopulationAllele::addDiseaseAlleles(const DBCitationMap& disease_allele_map) {

  for (auto const& [allele, pmid_set] : disease_allele_map) {

    auto publication_map = pubmed_requestor_ptr_->getCachedPublications(pmid_set);

    static const PlasmodiumFilter pf_filter;
    std::set<std::string> filtered_pmid;

    for (auto const& [pmid, publication_ptr] : publication_map) {

      if (pf_filter.applyFilter(*publication_ptr)) {

        filtered_pmid.insert(pmid);

      }

    }

    if (not filtered_pmid.empty()) {

      auto [iter, result] = disease_allele_map_.try_emplace(allele, filtered_pmid);
      if (not result) {

        ExecEnv::log().error("GeneratePopulationAllele::addDiseaseAlleles; unexpected duplicate allele: {}", allele);

      }

    }

  }

}


void kga::GeneratePopulationAllele::processPopulation( const std::shared_ptr<const PopulationDB>& population_ptr,
                                                       const std::shared_ptr<const SortedVariantAnalysis>& sorted_variant_ptr) {

  ExecEnv::log().info("Begin analyzing Literature Population: {}, with Genomes: {}", population_ptr->populationId(), population_ptr->getMap().size());

  WorkflowThreads thread_pool(WorkflowThreads::defaultThreads());
  std::vector<std::future<ThreadReturnType>> future_vector;

  // Create a disease allele map resource.
  std::shared_ptr<const DBCitationMap> disease_allele_ptr(std::make_shared<const DBCitationMap>(disease_allele_map_));

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    // Function ptr, args by value.
    std::future<std::pair<std::string, std::set<std::string>>> future = thread_pool.enqueueFuture(&GeneratePopulationAllele::getGenomePublications,
                                                                                                  genome_ptr,
                                                                                                  disease_allele_ptr);
    future_vector.push_back(std::move(future));

  }

  // Unpack the results.
  for (auto& future : future_vector) {

    auto [genome_id, allele_set] = future.get();

    reference_ethnic_.genomeAnalysis(genome_id, 1, genome_aux_ptr_);

    for (auto const& allele : allele_set) {

      auto result = variant_allele_map_.find(allele);
      if (result == variant_allele_map_.end()) {

        variant_allele_map_.emplace(allele, std::set<std::string>{genome_id});

      } else {

        auto& [rs_key, genome_set] = *result;
        genome_set.insert(genome_id);

      }

    }

  }

  // Setup the ensembl codes.
  for (auto const& [allele, genome_set] : variant_allele_map_) {

    auto result = sorted_variant_ptr->alleleEnsemblMap()->find(allele);
    if (result != sorted_variant_ptr->alleleEnsemblMap()->end()) {

      auto const& [rs_key, ensembl_codes] = *result;
      allele_ensembl_codes_.emplace(allele, ensembl_codes);

    }

  }

  ExecEnv::log().info("Completed analyzing Literature Population: {}, with Genomes: {}", population_ptr->populationId(), population_ptr->getMap().size());

}

kga::GeneratePopulationAllele::ThreadReturnType kga::GeneratePopulationAllele::getGenomePublications( std::shared_ptr<const GenomeDB> genome_ptr,
                                                                                                      std::shared_ptr<const DBCitationMap> disease_cited_alleles_ptr) {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Local class to process all the genome data.
  class ProcessGenome {

  public:

    explicit ProcessGenome(std::shared_ptr<const DBCitationMap> disease_cited_alleles_ptr) : disease_cited_alleles_ptr_(std::move(disease_cited_alleles_ptr)) {}
    ~ProcessGenome() = default;

    bool forEachVariant(const std::shared_ptr<const Variant>& variant_ptr) {

      if (disease_cited_alleles_ptr_->contains(variant_ptr->identifier())) {

        variant_allele_set_.insert(variant_ptr->identifier());

      }

      return true;

    }

    std::set<std::string> getAlleleSet() { return variant_allele_set_; }

  private:

    std::shared_ptr<const DBCitationMap> disease_cited_alleles_ptr_;
    std::set<std::string> variant_allele_set_;

  };
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


  ProcessGenome process_genome(std::move(disease_cited_alleles_ptr));

  genome_ptr->processAll(process_genome, &ProcessGenome::forEachVariant);

  return { genome_ptr->genomeId(), process_genome.getAlleleSet() };

}


void kga::GeneratePopulationAllele::writePopLiterature(const std::string& output_file) const {


  std::ofstream out_file(output_file);

  if (not out_file.good()) {

    ExecEnv::log().error("GeneratePopulationAllele::writePopLiterature; cannot open output file: {}", output_file);
    return;

  }

  ExecEnv::log().info("Writing literature analysis for: {} variants to file: {}", variant_allele_map_.size(), output_file);

  // Sort by number of genomes.
  std::multimap<size_t, std::string> allele_pubcount_map;
  for (auto const& [allele_rs_id, genome_set] : variant_allele_map_) {

    allele_pubcount_map.emplace(genome_set.size(), allele_rs_id);

  }

  for (auto iter = allele_pubcount_map.rbegin(); iter != allele_pubcount_map.rend(); ++iter) {

    // Unwrap the variables.
    auto const& [genome_count, rs_id] = *iter;

    auto result = variant_allele_map_.find(rs_id);
    if (result == variant_allele_map_.end()) {

      ExecEnv::log().error("GeneratePopulationAllele::writePopLiterature; allele id: {} not found in genome map", rs_id);
      continue;

    }

    // Calculate the ethnic breakdown.
    auto const& [rs_key, genome_set] = *result;
    GeneEthnicitySex allele_ethnic;
    allele_ethnic.updatePopulations(genome_aux_ptr_);
    for (auto const& genome : genome_set) {

      allele_ethnic.genomeAnalysis(genome, 1 , genome_aux_ptr_);

    }
    if (reference_ethnic_.superPopulation().size() != allele_ethnic.superPopulation().size()) {

      ExecEnv::log().error("GeneratePopulationAllele::writePopLiterature; allele id: {}, reference super populations: {}, allele super populations: {}",
                           rs_id, reference_ethnic_.superPopulation().size(), allele_ethnic.superPopulation().size());
      continue;

    }

    out_file << "\n******************************************\n\n";
    double pop_freq = (static_cast<double>(genome_count) * 100.0) / static_cast<double>(reference_ethnic_.total());
    out_file << "Genome Count: " << genome_count << '/'  << reference_ethnic_.total() << " (" << pop_freq << "%)\n\n";

    auto ref_iter = reference_ethnic_.superPopulation().begin();
    auto allele_iter = allele_ethnic.superPopulation().begin();

    while (ref_iter != reference_ethnic_.superPopulation().end() and allele_iter != allele_ethnic.superPopulation().end()) {

      auto const& [ref_pop, ref_count] = *ref_iter;
      auto const& [allele_pop, allele_count] = *allele_iter;
      double pop_percent = (static_cast<double>(allele_count) * 100.0) / static_cast<double>(ref_count);

      out_file << ref_pop << " " << allele_count << "/" << ref_count << " (";
      out_file << pop_percent << "%), ";

      ++ref_iter;
      ++allele_iter;

    }
    out_file << "\n\n";

    ref_iter = reference_ethnic_.population().begin();
    allele_iter = allele_ethnic.population().begin();

    while (ref_iter != reference_ethnic_.population().end() and allele_iter != allele_ethnic.population().end()) {

      auto const& [ref_pop, ref_count] = *ref_iter;
      auto const& [allele_pop, allele_count] = *allele_iter;
      double pop_percent = (static_cast<double>(allele_count) * 100.0) / static_cast<double>(ref_count);

      out_file << ref_pop << " " << allele_count << "/" << ref_count << " (";
      out_file << pop_percent << "%), ";

      ++ref_iter;
      ++allele_iter;

    }
    out_file << "\n\n";


    auto find_ensembl = allele_ensembl_codes_.find(rs_id);
    if (find_ensembl != allele_ensembl_codes_.end()) {

      auto const& [rs_id, ensembl_set] = *find_ensembl;
      auto [concat_symbol, concat_id] = generateGeneCodes(ensembl_set);
      out_file << rs_id << "|" << concat_symbol << "|" << concat_id << '\n';

    } else {

      out_file << rs_id << "|" << "***No Gene***" << "|" << "***No Gene***" << '\n';

    }

    out_file << "\n\n******************************************" << '\n';

    // print all the publications
    auto pmid_set = getDiseaseCitations(rs_id);
    std::vector<std::string> pmid_vector;
    for (auto const& pmid : pmid_set) {

      pmid_vector.push_back(pmid);

    }

    // Get the literature for this allele;
    auto literature_map = pubmed_requestor_ptr_->getCachedPublications(pmid_vector);

    for (auto const& [pmid, publication_ptr] : literature_map) {

      out_file << '\n';
      publication_ptr->extendedBiblio(out_file);
      out_file << '\n';

    }

  } // for all variants.

}


std::set<std::string> kga::GeneratePopulationAllele::getDiseaseCitations(const std::string& rs_identifier) const {

  std::set<std::string> allele_disease_set;

  auto result = disease_allele_map_.find(rs_identifier);
  if (result != disease_allele_map_.end()) {

    auto const& [rs_key, pmid_vector] = *result;
    for (auto const& pmid : pmid_vector) {

      allele_disease_set.insert(pmid);

    }

  }

  return allele_disease_set;

}


std::pair<std::string, std::string> kga::GeneratePopulationAllele::generateGeneCodes(const std::set<std::string>& gene_id_set) const {


  std::set<std::string> symbol_set;
  for (auto const& allele_id : gene_id_set) {

    auto symbol_vector = uniprot_nomenclature_ptr_->ensemblToSymbol(allele_id);
    if (not symbol_vector.empty()) {

      for (auto const& symbol : symbol_vector) {

        if (not symbol.empty()) {

          symbol_set.insert(symbol);

        }

      }

    }

  }

  // Some of the ids may be Entrez.
  for (auto const& allele_id : gene_id_set) {

    std::string symbol = entrez_nomenclature_ptr_->entrezToSymbol(allele_id);
    if (not symbol.empty()) {

      symbol_set.insert(symbol);

    }

  }

  // Generate Id and symbol strings
  std::string concat_id;
  for (auto const& allele_id : gene_id_set) {

    concat_id += allele_id;
    if (allele_id != *gene_id_set.rbegin()) {

      concat_id += CONCAT_DELIMITER_;

    }

  }

  std::string concat_symbol;
  for (auto const& symbol : symbol_set) {

    concat_symbol += symbol;
    if (symbol != *symbol_set.rbegin()) {

      concat_symbol += CONCAT_DELIMITER_;

    }

  }

  return {concat_symbol, concat_id};

}