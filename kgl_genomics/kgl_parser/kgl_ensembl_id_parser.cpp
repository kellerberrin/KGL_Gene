//
// Created by kellerberrin on 16/3/21.
//

#include "kgl_ensembl_id_parser.h"
#include "kel_utility.h"



namespace kgl = kellerberrin::genome;


void kgl::EnsemblHGNCResource::IndexHGNC() {

  static bool error{false};

  // First load up the canonical records.
  for (auto const& [hgnc_id, canonical, ensembl_id] : synonym_vector_) {

    if (hgnc_id.empty() or canonical.empty()) {

      continue;

    }
    auto [iter, result] = hgnc_emsembl_map_.try_emplace(hgnc_id, ensembl_id);
    if (not result and not error) {

      ExecEnv::log().warn("EnsemblHGNCResource::IndexHGNC; duplicate HGNC records ({})", hgnc_id);
      error = true;

    }

  }
  // Then load up all other records.
  for (auto const& [hgnc_id, canonical, ensembl_id] : synonym_vector_) {

    if (hgnc_id.empty() or not canonical.empty()) {

      continue;

    }

    hgnc_emsembl_map_.emplace(hgnc_id, ensembl_id);

  }

  ExecEnv::log().info("EnsemblHGNCResource loaded {}, HGNC, Ensembl Id lookup pairs", hgnc_emsembl_map_.size());

}


[[nodiscard]] std::string kgl::EnsemblHGNCResource::HGNCToEnsembl(const std::string& hgnc_id) const {

  auto result = hgnc_emsembl_map_.find(hgnc_id);
  if (result == hgnc_emsembl_map_.end()) {

    return "";

  }

  auto [hgnc, ensembl] = *result;

  return ensembl;

}


[[nodiscard]] bool kgl::ParseGeneIdents::parseIdentFile(const std::string& file_name) {


  auto parsed_record_ptr = parseFlatFile(file_name, SquareTextParser::DELIMITER_CSV);

  if (parsed_record_ptr->getRowVector().size() < MINIMUM_ROW_COUNT_) {

    ExecEnv::log().error("ParseGeneIdents::parseUniprotFile; Row count: {} for file: {} is below minimum",
                         parsed_record_ptr->getRowVector().size(), file_name);
    return false;

  }

  if (not parsed_record_ptr->checkRowSize(COLUMN_COUNT_)) {

    ExecEnv::log().error("ParseGeneIdents::parseUniprotFile; Not all rows have expected column count: {} for file: {}",
                         COLUMN_COUNT_, file_name);
    return false;

  }

  bool header_row{true};
  for (const auto& row_vector :  parsed_record_ptr->getRowVector()) {

    // Skip the header line.
    if (header_row) {

      header_row = false;
      continue;

    }

    GeneIDSynonyms gene_id_synonyms;
    gene_id_synonyms.HGNC_id = Utility::trimEndWhiteSpace(row_vector[HGNC_OFFSET_]);
    gene_id_synonyms.Canonical = Utility::trimEndWhiteSpace(row_vector[CANONICAL_OFFSET_]);
    gene_id_synonyms.ensembl_id = Utility::trimEndWhiteSpace(row_vector[ENSEMBL_OFFSET_]);
    synonym_vector_.push_back(gene_id_synonyms);

  }

  ExecEnv::log().info("ParseGeneIdents::parseUniprotFile; Parsed: {} HGNC and Ensembl gene identifiers from file: {}",
                      synonym_vector_.size() , file_name);
  return true;

}
