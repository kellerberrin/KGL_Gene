//
// Created by kellerberrin on 2/8/21.
//

#include "kgl_entrez_parser.h"
#include "kel_utility.h"



namespace kgl = kellerberrin::genome;


void kgl::EntrezResource::IndexEntrez() {

  static bool error{false};

  // First load up the canonical records.
  for (auto const& entrez_record : entrez_vector_) {

    if (entrez_record.symbol_id.empty()) {

      continue;

    }
    auto [iter, result] = entrez_map_.try_emplace(entrez_record.symbol_id, entrez_record);
    if (not result and not error) {

      ExecEnv::log().warn("EntrezResource::IndexEntrez; duplicate Symbol record ({})",  entrez_record.symbol_id);
      error = true;

    }

  }

  ExecEnv::log().info("EntrezResource loaded {}, (Symbol, Entrez) Id lookup pairs", entrez_map_.size());

}


[[nodiscard]] std::string kgl::EntrezResource::symbolToEntrez(const std::string& symbol_id) const {

  auto result = entrez_map_.find(symbol_id);
  if (result == entrez_map_.end()) {

    return "";

  }

  auto [symbol, entrez_record] = *result;

  return entrez_record.entrez_id;

}


[[nodiscard]] bool kgl::ParseEntrez::parseEntrezFile(const std::string& file_name) {


  auto parsed_record_ptr = parseFlatFile(file_name, SquareTextParser::DELIMITER_TSV);

  if (parsed_record_ptr->getRowVector().size() < MINIMUM_ROW_COUNT_) {

    ExecEnv::log().error("ParseEntrez::parseEntrezFile; Row count: {} for file: {} is below minimum",
                         parsed_record_ptr->getRowVector().size(), file_name);
    return false;

  }

  if (not parsed_record_ptr->checkRowSize(COLUMN_COUNT_)) {

    ExecEnv::log().error("ParseEntrez::parseEntrezFile; Not all rows have expected column count: {} for file: {}",
                         COLUMN_COUNT_, file_name);
    return false;

  }

  ExecEnv::log().info("Begin Parsing Entrez Gene Resource for file: {}", file_name);

  // Header line is prefixed with '#' and is automatically stripped off.
  for (const auto& row_vector :  parsed_record_ptr->getRowVector()) {

    EntrezRecord entrez_record;
    entrez_record.entrez_id = Utility::trimEndWhiteSpace(row_vector[ENTREZ_OFFSET_]);
    entrez_record.symbol_id = Utility::trimEndWhiteSpace(row_vector[SYMBOL_OFFSET_]);
    entrez_record.description = Utility::trimEndWhiteSpace(row_vector[DESCRIPTION_OFFSET_]);
    entrez_vector_.push_back(entrez_record);

  }

  ExecEnv::log().info("ParseEntrez::parseEntrezFile; Parsed: {} Symbol, Entrez gene identifiers from file: {}",
                      entrez_vector_.size() , file_name);
  return true;

}
