//
// Created by kellerberrin on 31/8/21.
//

#include "kgl_analysis_mutation_gene_allele_pop.h"


namespace kgl = kellerberrin::genome;


void kgl::GeneratePopulationAllele::initialize( const std::shared_ptr<const HsGenomeAux>& genome_aux_ptr,
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

}


void kgl::GeneratePopulationAllele::processPopulation(const std::shared_ptr<const PopulationDB>& population_ptr) {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Local class to process all the population data.
  class ProcessPopulation {

  public:

    explicit ProcessPopulation(DBCitationMap disease_allele_map) : disease_allele_map_(std::move(disease_allele_map)) {}
    ~ProcessPopulation() = default;

    bool forEachVariant(const std::shared_ptr<const Variant>& variant_ptr) {

      if (disease_allele_map_.contains(variant_ptr->identifier())) {

        auto result = variant_allele_map_.find(variant_ptr->identifier());
        if (result == variant_allele_map_.end()) {

          variant_allele_map_.emplace(variant_ptr->identifier(), 1);

        } else {

          auto& [rs_key, genome_count] = *result;
          ++genome_count;

        }

      }

      return true;

    }

    VariantCountMap getCountMap() { return variant_allele_map_; }

  private:

    DBCitationMap disease_allele_map_;
    std::map<std::string, size_t> variant_allele_map_;

  };
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ExecEnv::log().info("Begin analyzing Literature Population: {}, with Genomes: {}", population_ptr->populationId(), population_ptr->getMap().size());
  // Check each variant.
  ProcessPopulation process_population(disease_allele_map_);

  population_ptr->processAll(process_population, &ProcessPopulation::forEachVariant);

  variant_allele_map_ = process_population.getCountMap();

  ExecEnv::log().info("Completed analyzing Literature Population: {}, with Genomes: {}", population_ptr->populationId(), population_ptr->getMap().size());

}


void kgl::GeneratePopulationAllele::processPopulationMT(const std::shared_ptr<const PopulationDB>& population_ptr) {

  ExecEnv::log().info("Begin analyzing Literature Population: {}, with Genomes: {}", population_ptr->populationId(), population_ptr->getMap().size());

  ThreadPool thread_pool(ThreadPool::defaultThreads());
  std::vector<std::future<ThreadReturnType>> future_vector;

  // Create a disease allele map resource.
  std::shared_ptr<const DBCitationMap> disease_allele_ptr(std::make_shared<const DBCitationMap>(disease_allele_map_));

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    // function, object_ptr, arg1
    std::future<std::pair<std::string, std::set<std::string>>> future = thread_pool.enqueueTask( &GeneratePopulationAllele::getGenomePublications,
                                                                                                  genome_ptr,
                                                                                                  disease_allele_ptr);
    future_vector.push_back(std::move(future));

  }

  // Unpack the results.
  variant_allele_map_.clear();
  for (auto& future : future_vector) {

    auto [genome_id, allele_set] = future.get();

    for (auto const& allele : allele_set) {

      auto result = variant_allele_map_.find(allele);
      if (result == variant_allele_map_.end()) {

        variant_allele_map_.emplace(allele, 1);

      } else {

        auto& [rs_key, count] = *result;
        ++count;

      }

    }

  }

  ExecEnv::log().info("Completed analyzing Literature Population: {}, with Genomes: {}", population_ptr->populationId(), population_ptr->getMap().size());

}

kgl::GeneratePopulationAllele::ThreadReturnType kgl::GeneratePopulationAllele::getGenomePublications( std::shared_ptr<const GenomeDB> genome_ptr,
                                                                                                      std::shared_ptr<const DBCitationMap> disease_cited_alleles) {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Local class to process all the genome data.
  class ProcessGenome {

  public:

    explicit ProcessGenome(std::shared_ptr<const DBCitationMap> disease_cited_alleles) : disease_cited_alleles_(std::move(disease_cited_alleles)) {}
    ~ProcessGenome() = default;

    bool forEachVariant(const std::shared_ptr<const Variant>& variant_ptr) {

      if (disease_cited_alleles_->contains(variant_ptr->identifier())) {

        variant_allele_set_.insert(variant_ptr->identifier());

      }

      return true;

    }

    std::set<std::string> getAlleleSet() { return variant_allele_set_; }

  private:

    std::shared_ptr<const DBCitationMap> disease_cited_alleles_;
    std::set<std::string> variant_allele_set_;

  };
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


  ProcessGenome process_genome(std::move(disease_cited_alleles));

  genome_ptr->processAll(process_genome, &ProcessGenome::forEachVariant);

  return { genome_ptr->genomeId(), process_genome.getAlleleSet() };

}


void kgl::GeneratePopulationAllele::writeOutput(const std::string& output_file) const {


  std::ofstream out_file(output_file);

  if (not out_file.good()) {

    ExecEnv::log().error("GeneratePopulationAllele::writeOutput; cannot open output file: {}", output_file);
    return;

  }

  ExecEnv::log().info("Writing literature analysis for: {} variants to file: {}", variant_allele_map_.size(), output_file);

  // Sort by number of genomes.
  std::multimap<size_t, std::string> allele_pubcount_map;
  for (auto const& [allele_rs_id, genome_count] : variant_allele_map_) {

    allele_pubcount_map.emplace(genome_count, allele_rs_id);

  }

  for (auto iter = allele_pubcount_map.rbegin(); iter != allele_pubcount_map.rend(); ++iter) {

    // Unwrap the variables.
    auto const& [genome_count, rs_id] = *iter;
    std::vector<std::string> allele_id_array{rs_id};
    auto [concat_symbol, concat_id] = generateGeneCodes(allele_id_array);

    out_file << "\n******************************************\n\n";
    out_file << "Genome Count: " << genome_count << '\n';
    out_file << rs_id << "|" << concat_symbol << "|" << concat_id << '\n';
    out_file << "\n\n******************************************" << '\n';

    // print all the publications
    auto pmid_set = getDiseaseCitations(rs_id);
    std::vector<std::string> pmid_vector;
    for (auto const& pmid : pmid_set) {

      pmid_vector.push_back(pmid);

    }

    // Get the literature for this allele;
    auto literature_map = pubmed_requestor_ptr_->getCachedPublications(pmid_vector);

    for (auto const& [pmid, publication] : literature_map) {

      out_file << '\n';
      publication.extendedBiblio(out_file);
      out_file << '\n';

    }

  } // for all variants.

}


std::set<std::string> kgl::GeneratePopulationAllele::getCitations(const std::string& rs_identifier) const {

  std::set<std::string> allele_pmid_set;
  if (not rs_identifier.empty()) {

    auto find_result = disease_allele_map_.find(rs_identifier);
    if (find_result != disease_allele_map_.end()) {

      auto const& [rsid, citations] = *find_result;
      allele_pmid_set = citations;

    }

  }

  return allele_pmid_set;

}


bool kgl::GeneratePopulationAllele::citationsExist(const std::string& rs_identifier) const {

  if (not rs_identifier.empty()) {

    return allele_citation_ptr_->alleleIndexedCitations().contains(rs_identifier);

  }

  return false;

}


std::set<std::string> kgl::GeneratePopulationAllele::getDiseaseCitations(const std::string& rs_identifier) const {

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


std::pair<std::string, std::string> kgl::GeneratePopulationAllele::generateGeneCodes(const std::vector<std::string>& ensembl_entrez_codes) const {


  std::set<std::string> gene_id_set;
  for (auto const& allele_id : ensembl_entrez_codes) {

    if (not allele_id.empty()) {

      gene_id_set.insert(allele_id);

    }

  }

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

      concat_id += CONCATENATE_VEP_FIELDS_;

    }

  }

  std::string concat_symbol;
  for (auto const& symbol : symbol_set) {

    concat_symbol += symbol;
    if (symbol != *symbol_set.rbegin()) {

      concat_symbol += CONCATENATE_VEP_FIELDS_;

    }

  }

  return {concat_symbol, concat_id};

}