//
// Created by kellerberrin on 20/6/21.
//

#include "kgl_analysis_mutation_gene_allele.h"
#include "kgl_variant_db_freq.h"
#include "kgl_variant_sort.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kel_distribution.h"

#include <fstream>

namespace kgl = kellerberrin::genome;


void kgl::GenerateGeneAllele::initialize( const std::vector<std::string>& symbol_gene_list,
                                          const std::shared_ptr<const UniprotResource>& uniprot_nomenclature_ptr,
                                          const std::shared_ptr<const EntrezResource>& entrez_nomenclature_ptr,
                                          const std::shared_ptr<const CitationResource>& allele_citation_ptr,
                                          const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr) {


  uniprot_nomenclature_ptr_ = uniprot_nomenclature_ptr;
  entrez_nomenclature_ptr_ = entrez_nomenclature_ptr;
  allele_citation_ptr_ = allele_citation_ptr;
  pubmed_requestor_ptr_ = pubmed_requestor_ptr;

  for (auto const& symbol : symbol_gene_list) {

    auto ensembl_symbol_list = uniprot_nomenclature_ptr_->symbolToEnsembl(symbol);
    for (auto const& ensembl :  ensembl_symbol_list) {

      auto [iter, insert_result] = ensembl_symbol_map_.try_emplace(ensembl, symbol);
      if (not insert_result) {

        ExecEnv::log().error("GenerateGeneAllele::initialize; could not add (duplicate) gene ensembl code: {} for gene symbol: {}", ensembl, symbol);

      }

    }

  }

  cited_allele_map_.clear();

}


void kgl::GenerateGeneAllele::addGeneCitedVariants(const std::shared_ptr<const SortedVariantAnalysis>& sorted_variants) {

  std::vector<std::string> ensembl_list;
  for (auto const& [ensembl, symbol] : ensembl_symbol_map_) {

    ensembl_list.push_back(ensembl);

  }

  auto filtered_variants = sorted_variants->filterEnsembl(ensembl_list);

  for (auto const& [ensembl_id, variant_ptr] : filtered_variants) {

    if (not getCitations(variant_ptr->identifier()).empty()) {

      auto result = cited_allele_map_.find(variant_ptr->identifier());
      if (result == cited_allele_map_.end()) {

        std::pair<std::shared_ptr<const Variant>,std::vector<std::string>> value_pair{variant_ptr, std::vector<std::string>{ensembl_id}};
        cited_allele_map_.emplace(variant_ptr->identifier(), value_pair);

      } else {

        auto& [allele_rs_key, variant_code_pair] = *result;
        auto& [value_variant_ptr, id_array] = variant_code_pair;

        id_array.push_back(ensembl_id);

      }

    }

  }

}

void kgl::GenerateGeneAllele::addDiseaseCitedVariants(const std::shared_ptr<const SortedVariantAnalysis>& sorted_variants) {

  // To save space only add citations that are relevant to the disease MeSH code.
  for (auto const& [ensembl_id, variant_ptr] : *(sorted_variants->ensemblMap())) {

    if (not getDiseaseCitations(variant_ptr->identifier()).empty()) {

      auto result = cited_allele_map_.find(variant_ptr->identifier());
      if (result == cited_allele_map_.end()) {

        std::pair<std::shared_ptr<const Variant>,std::vector<std::string>> value_pair{variant_ptr, std::vector<std::string>{ensembl_id}};
        cited_allele_map_.emplace(variant_ptr->identifier(), value_pair);

      } else {

        auto& [allele_rs_key, variant_code_pair] = *result;
        auto& [value_variant_ptr, id_array] = variant_code_pair;

        id_array.push_back(ensembl_id);

      }

    }

  }

}


std::set<std::string> kgl::GenerateGeneAllele::getCitations(const std::string& rs_identifier) const {

  std::set<std::string> allele_pmid_set;
  if (not rs_identifier.empty()) {

    auto find_result = allele_citation_ptr_->citationMap().find(rs_identifier);
    if (find_result != allele_citation_ptr_->citationMap().end()) {

      auto const& [rsid, citations] = *find_result;
      std::string cite_string;
      for (auto const& cite : citations) {

        allele_pmid_set.insert(cite);

      }

    }

  }

  return allele_pmid_set;

}

std::set<std::string> kgl::GenerateGeneAllele::getDiseaseCitations(const std::string& rs_identifier) const {

  std::set<std::string> allele_disease_set;
  auto rs_pmid_set = getCitations(rs_identifier);

  for (auto const& pmid : rs_pmid_set) {

    if (disease_citations_.contains(pmid)) {

      allele_disease_set.insert(pmid);

    }

  }

  return allele_disease_set;

}


void kgl::GenerateGeneAllele::writeHeader(std::ofstream& outfile, char delimiter) {

  outfile << "SymbolGeneId" << delimiter
          << "EnsemblGeneId" << delimiter
          << "VariantId" << delimiter
          << "ContigId" << delimiter
          << "Offset" << delimiter
          << "Reference" << delimiter
          << "Alternate" << delimiter
          << "CiteCount" << delimiter
          << "DiseaseCites" << delimiter
          << "GlobalFreq" << delimiter
          << "AFRFreq"  << delimiter;

  for (auto const& field_id : VEP_FIELD_LIST_) {

    outfile << field_id << delimiter;

  }

  outfile << "TotalCount" << delimiter
          << "AltCount" << delimiter
          << "AFRCount" << delimiter
          << "AFRAltCount" << delimiter
          << "AFR_Freq%" << delimiter
          << "NonAFR_Freq%" << delimiter
          << "DiseasePMIDs" << '\n';

}


void kgl::GenerateGeneAllele::writeOutput(const std::string& output_file, char delimiter) const {

  std::ofstream out_file(output_file);

  if (not out_file.good()) {

    ExecEnv::log().error("GenerateGeneAllele::writeOutput; cannot open output file: {}", output_file);
    return;

  }

  ExecEnv::log().info("Writing VEP, PMID information for {} variants to file: {}", cited_allele_map_.size(), output_file);

  writeHeader(out_file, delimiter);

  // The key is the allele unique 'rs' code, the value pair is .first is a pointer to the variant, .second is a vector of all the gene codes for the allele.
  for (auto const& [allele_rs_id, variant_ptr_id_pair] : cited_allele_map_) {

    auto const& [variant_ptr, allele_id_array] = variant_ptr_id_pair;

    auto [concat_symbol, concat_id] = generateGeneCodes(allele_id_array);

    out_file << concat_symbol << delimiter;

    out_file << concat_id << delimiter
             << variant_ptr->identifier() << delimiter
             << variant_ptr->contigId() << delimiter
             << variant_ptr->offset() << delimiter
             << variant_ptr->reference().getSequenceAsString() << delimiter
             << variant_ptr->alternate().getSequenceAsString() << delimiter;

    auto pmid_set = getCitations(variant_ptr->identifier());
    out_file << pmid_set.size() << delimiter;

    auto disease_pmid_set = getDiseaseCitations(variant_ptr->identifier());
    out_file << disease_pmid_set.size() << delimiter;

    double global_freq{0.0};
    auto frequency_opt = FrequencyDatabaseRead::superPopFrequency(*variant_ptr, ALL_SUPER_POP_);
    if (frequency_opt) {

      global_freq = frequency_opt.value();

    }

    out_file << global_freq << delimiter;

    double afr_freq{0.0};
    auto afr_frequency_opt = FrequencyDatabaseRead::superPopFrequency(*variant_ptr, AFR_SUPER_POP_);
    if (afr_frequency_opt) {

      afr_freq = afr_frequency_opt.value();

    }

    out_file << afr_freq << delimiter;

    auto vep_fields =  retrieveVepFields( variant_ptr, VEP_FIELD_LIST_);
    for (auto const& [field_id, field_value] : vep_fields) {

      out_file << field_value << delimiter;

    }


    auto total_allele_count_opt =  FrequencyDatabaseRead::superPopTotalAlleles(*variant_ptr, ALL_SUPER_POP_);
    auto alt_allele_count_opt =  FrequencyDatabaseRead::superPopAltAlleles(*variant_ptr, ALL_SUPER_POP_);
    auto afr_total_allele_count_opt =  FrequencyDatabaseRead::superPopTotalAlleles(*variant_ptr, AFR_SUPER_POP_);
    auto afr_alt_allele_count_opt =  FrequencyDatabaseRead::superPopAltAlleles(*variant_ptr, AFR_SUPER_POP_);

    size_t total_allele_count{0};
    size_t alt_allele_count{0};
    size_t afr_total_allele_count{0};
    size_t afr_alt_allele_count{0};

    if (total_allele_count_opt and alt_allele_count_opt and afr_total_allele_count_opt and afr_alt_allele_count_opt) {

      total_allele_count = total_allele_count_opt.value();
      alt_allele_count = alt_allele_count_opt.value();
      afr_total_allele_count = afr_total_allele_count_opt.value();
      afr_alt_allele_count = afr_alt_allele_count_opt.value();

    }

    out_file << total_allele_count << delimiter
             << alt_allele_count << delimiter
             << afr_total_allele_count << delimiter
             << afr_alt_allele_count << delimiter;

    if (afr_total_allele_count > 0) {

      double afr_percent = (static_cast<double>(afr_alt_allele_count) / static_cast<double>(afr_total_allele_count)) * 100.0;
      out_file << afr_percent << delimiter;

    } else {

      out_file << 0 << delimiter;

    }

    int64_t non_afr_count = total_allele_count - afr_total_allele_count;
    if (non_afr_count  > 0) {

      double non_afr_percent = ((static_cast<double>(alt_allele_count)-static_cast<double>(afr_alt_allele_count)) / static_cast<double>(non_afr_count)) * 100.0;
      out_file << non_afr_percent  << delimiter;

    } else {

      out_file << 0 << delimiter;

    }

    std::string disease_cite_string;
    for (auto const& pmid : disease_pmid_set) {

      disease_cite_string += pmid;
      if (pmid != *pmid_set.rbegin()) {

        disease_cite_string += CONCATENATE_VEP_FIELDS_;

      }

    }
    out_file << disease_cite_string << '\n';

  }

}


std::map<std::string, std::string> kgl::GenerateGeneAllele::retrieveVepFields( const std::shared_ptr<const Variant>& variant_ptr,
                                                                               const std::vector<std::string>& field_list) {

  static std::vector<std::pair<std::string, size_t>> field_indicies;
  static bool initialized{false};
  static bool error_flag{false};

  if (not initialized) {

    field_indicies = InfoEvidenceAnalysis::getVepIndexes(*variant_ptr, field_list);
    initialized = true;

  }

  auto field_vector = InfoEvidenceAnalysis::getVepData(*variant_ptr, field_indicies);

  std::map<std::string, std::set<std::string>> aggregated_fields;
  // There may be multiple VEP fields.
  for (auto& field : field_vector) {

    if (field.size() != field_list.size()) {

      if (not error_flag) {

        ExecEnv::log().warn("GenerateGeneAllele::retrieveVepFields; requested fields: {}, returned VEP fields: {}", field_list.size(), field.size());
        error_flag = true;

      }

      continue;

    }

    for (auto const& [field_id, field_value] : field) {

      auto result = aggregated_fields.find(field_id);
      if (result == aggregated_fields.end()) {

        aggregated_fields.emplace(field_id, std::set<std::string>{field_value});

      } else {

        if (not field_value.empty()) {

          auto& [agg_field_id, agg_set_value] = *result;
          agg_set_value.insert(field_value);

        }

      }

    }

  }

  // Concatenate multiple field values.
  std::map<std::string, std::string> vep_fields;
  for (auto const&[field_id, field_value_set] : aggregated_fields) {

    std::string concat_fields;
    for (auto const& field_value : field_value_set) {

      concat_fields += field_value + CONCATENATE_VEP_FIELDS_;

    }

    vep_fields[field_id] = concat_fields;

  }

  return vep_fields;

}


void kgl::GenerateGeneAllele::writeLiteratureAlleleSummary(const std::string& output_file) {

  std::ofstream out_file(output_file);

  if (not out_file.good()) {

    ExecEnv::log().error("GenerateGeneAllele::writeLiteratureAlleleSummary; cannot open output file: {}", output_file);
    return;

  }

  ExecEnv::log().info("Writing literature  summaries for: {} variants to file: {}", cited_allele_map_.size(), output_file);

  std::vector<std::string> pmid_vector;
  std::map<std::string, std::vector<std::pair<std::shared_ptr<const Variant>, std::vector<std::string>>>> pmid_variant_map;
  for (auto const& [allele_rs_id, variant_ptr_id_pair] : cited_allele_map_) {

    auto const& [variant_ptr, allele_id_array] = variant_ptr_id_pair;

    auto pmid_set = getDiseaseCitations(variant_ptr->identifier());
    for (auto const& pmid : pmid_set) {

      pmid_vector.push_back(pmid);
      auto result = pmid_variant_map.find(pmid);
      if (result == pmid_variant_map.end()) {

        std::vector<std::pair<std::shared_ptr<const Variant>, std::vector<std::string>>> variant_vector;
        variant_vector.emplace_back(variant_ptr, allele_id_array);
        pmid_variant_map.emplace(pmid, variant_vector);

      } else {

        auto& [pmid_key, variant_vector] = *result;
        variant_vector.emplace_back(variant_ptr, allele_id_array);

      }

    }

  }

  auto literature_map = pubmed_requestor_ptr_->getCachedPublications(pmid_vector);

  // Resort the literature map by number of citations.
  std::multimap<size_t, PublicationSummary> citation_rank_map;
  for (auto const& [pmid, publication] : literature_map) {

    citation_rank_map.emplace(publication.citedBy().size(), publication);

  }

  // Print most cited articles first.
  for (auto iter = citation_rank_map.rbegin(); iter != citation_rank_map.rend(); ++iter) {

    auto const& [cites, publication] = *iter;

    out_file << "\n******************************************" << '\n';
    auto result = pmid_variant_map.find(publication.pmid());
    if (result == pmid_variant_map.end()) {

      ExecEnv::log().error("GenerateGeneAllele::writeLiteratureAlleleSummary; no alleles found for publication pmid: {}", publication.pmid());

    } else {

      auto const& [pmid_key, variant_vector] = *result;

      for (auto const& [variant_ptr, allele_id_array] : variant_vector) {

        auto [concat_symbol, concat_id] = generateGeneCodes(allele_id_array);

        out_file << variant_ptr->identifier() << "|" << concat_symbol << "|" << concat_id << '\n';

      } // For all alleles referenced by this publication.

      out_file << "******************************************" << '\n';

      out_file << '\n';
      publication.extendedBiblio(out_file);
      out_file << '\n';

    } // if found publication.

  } // citation ranked publication.

}


// For each Allele print all the relevant publications.
void kgl::GenerateGeneAllele::writeAlleleLiteratureSummary(const std::string& output_file) {

  std::ofstream out_file(output_file);

  if (not out_file.good()) {

    ExecEnv::log().error("GenerateGeneAllele::writeAlleleLiteratureSummary; cannot open output file: {}", output_file);
    return;

  }

  ExecEnv::log().info("Writing literature  summaries for: {} variants to file: {}", cited_allele_map_.size(), output_file);

  // Sort by number of publications.
  std::multimap<size_t, std::pair<std::shared_ptr<const Variant>, std::vector<std::string>>> allele_pubcount_map;
  for (auto const& [allele_rs_id, variant_ptr_id_pair] : cited_allele_map_) {

    auto const& [variant_ptr, allele_id_array] = variant_ptr_id_pair;

    auto pmid_set = getDiseaseCitations(variant_ptr->identifier());

    std::pair<std::shared_ptr<const Variant>, std::vector<std::string>> value_pair{variant_ptr, allele_id_array};
    allele_pubcount_map.emplace(pmid_set.size(), value_pair);

  }

  for (auto iter = allele_pubcount_map.rbegin(); iter != allele_pubcount_map.rend(); ++iter) {

    // Unwrap the variables.
    auto const& [pub_count, variant_ptr_id_pair] = *iter;
    auto const& [variant_ptr, allele_id_array] = variant_ptr_id_pair;

    auto [concat_symbol, concat_id] = generateGeneCodes(allele_id_array);

    out_file << "\n******************************************\n\n";
    out_file << variant_ptr->identifier() << "|" << concat_symbol << "|" << concat_id << '\n';
    out_file << "\n\n******************************************" << '\n';

    // print all the publications
    auto pmid_set = getDiseaseCitations(variant_ptr->identifier());
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



std::pair<std::string, std::string> kgl::GenerateGeneAllele::generateGeneCodes(const std::vector<std::string>& ensembl_entrez_codes) const {



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