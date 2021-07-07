//
// Created by kellerberrin on 2/7/21.
//

#include "kel_utility.h"
#include "kgl_uniprot_parser.h"



namespace kgl = kellerberrin::genome;


void kgl::UniprotResource::createIndexes() {

  size_t info_count{0};
  for (auto const& [uniprot_id, info_map] : uniprot_map_) {

    uniprot_index_.emplace(info_map.uniprotac_id, info_map.uniprotkb_id);
    info_count += info_map.attribute_map.size();

    for (auto const& [field_id, field_value] : info_map.attribute_map) {

      if (field_id == GENE_NAME ) {

        symbol_index_.emplace(field_value, uniprot_id);

      } else if (field_id == GENE_SYNONYM) {

        synonym_index_.emplace(field_value, uniprot_id);

      }  else if (field_id == ENSEMBL_FIELD) {

        ensembl_index_.emplace(field_value, uniprot_id);

      }  else if (field_id == HGNC_FIELD) {

        hgnc_index_.emplace(field_value, uniprot_id);

      }  else if (field_id == ENTREZ_GENE) {

        entrez_index_.emplace(field_value, uniprot_id);

      }

    }

  } // for uniprot id.

  ExecEnv::log().info("Uniprot Resources Indexed; UniprotKB Id Count: {}, Info Records: {}, UniprotAC: {}, Symbol: {}, Synonym: {}, Ensembl: {}, HGNC: {}, Entrez: {}",
                       uniprot_map_.size(), info_count, uniprot_index_.size(), symbol_index_.size(), synonym_index_.size(),
                       ensembl_index_.size(), hgnc_index_.size(), entrez_index_.size());


}


std::vector<std::string> kgl::UniprotResource::lookupInfo( const std::string& lookup_value,
                                                           const std::multimap<std::string, std::string>& lookup_map,
                                                           const std::string& lookup_type) const {


  auto lower_lookup = lookup_map.lower_bound(lookup_value);
  auto const upper_lookup = lookup_map.upper_bound(lookup_value);

  std::vector<std::string> field_values;

  while(lower_lookup != upper_lookup) {

    auto const& [lookup, uniprotkb_id] = *lower_lookup;

    auto result = uniprot_map_.find(uniprotkb_id);
    if (result != uniprot_map_.end()) {

      auto const& [uniprot_key, attribute_map] = *result;
      auto lower_field = attribute_map.attribute_map.lower_bound(lookup_type);
      auto const upper_field = attribute_map.attribute_map.upper_bound(lookup_type);
      while(lower_field != upper_field) {

        auto const& [field_id, field_value] = *lower_field;
        field_values.push_back(field_value);
        ++lower_field;

      }

    }

    ++lower_lookup;

  }

  return field_values;

}

std::vector<std::string> kgl::UniprotResource::symbolToUniprot(const std::string& symbol) const {

  std::vector<std::string> uniprot;
  auto lower_lookup = symbol_index_.lower_bound(symbol);
  auto const upper_lookup = symbol_index_.upper_bound(symbol);

  while(lower_lookup != upper_lookup) {

    auto const&[sym_key, uniprotkb_id] = *lower_lookup;

    auto result = uniprot_map_.find(uniprotkb_id);
    if (result != uniprot_map_.end()) {

      auto const& [uniprot_key, attribute_map] = *result;
      uniprot.push_back(attribute_map.uniprotac_id);

    }

    ++lower_lookup;

  }

  if (uniprot.empty()) {

    uniprot = lookupInfo( symbol, synonym_index_, GENE_SYNONYM);

  }

  return uniprot;

}



std::vector<std::string> kgl::UniprotResource::HGNCToUniprot(const std::string& symbol) const {

  std::vector<std::string> uniprot;
  auto lower_lookup = hgnc_index_.lower_bound(symbol);
  auto const upper_lookup = hgnc_index_.upper_bound(symbol);

  while(lower_lookup != upper_lookup) {

    auto const&[sym_key, uniprotkb_id] = *lower_lookup;

    auto result = uniprot_map_.find(uniprotkb_id);
    if (result != uniprot_map_.end()) {

      auto const& [uniprot_key, attribute_map] = *result;
      uniprot.push_back(attribute_map.uniprotac_id);

    }

    ++lower_lookup;

  }

  return uniprot;

}


[[nodiscard]] bool kgl::ParseUniprotId::parseUniprotFile(const std::string& file_name) {


  auto parsed_record_ptr = parseFlatFile(file_name, SquareTextParser::DELIMITER_TSV);


  if (not parsed_record_ptr->checkRowSize(COLUMN_COUNT_)) {

    ExecEnv::log().error("ParseUniprotId::parseUniprotFile; Not all rows have expected column count: {} for file: {}",
                         COLUMN_COUNT_, file_name);
    return false;

  }

  // Simple state parser.
  enum class PARSER_STATE { INITIAL, PROCESS_FIELD};
  PARSER_STATE parser_state = PARSER_STATE::INITIAL;

  UniprotAttributeMap attribute_map;
  for (const auto& row_vector :  parsed_record_ptr->getRowVector()) {

    if (row_vector[FIELD_OFFSET_] == UniprotResource::UNIPROTKB_ID) {

      if (parser_state == PARSER_STATE::PROCESS_FIELD) {

        if (not attribute_map.uniprotkb_id.empty() and not attribute_map.attribute_map.empty()) {

          UniprotAttributeMap insert_map;
          insert_map.uniprotkb_id = attribute_map.uniprotkb_id;
          insert_map.uniprotac_id = attribute_map.uniprotac_id;

          auto [insert_iter, result] = uniprot_info_.try_emplace(attribute_map.uniprotkb_id, insert_map);
          if (not result) {

            ExecEnv::log().error("ParseUniprotId::parseUniprotFile, Unable to insert (duplicate) uniprotkb id: {}", attribute_map.uniprotkb_id);

          } else {

            // Store the attribute map efficiently
            auto& [uniprot_id, attrib_values] = *insert_iter;
            attrib_values.attribute_map = std::move(attribute_map.attribute_map);

          }

        } else {

          ExecEnv::log().error("ParseUniprotId::parseUniprotFile, Unexpected empty uniprotkb id, {field, value} record count: {}",
                               attribute_map.attribute_map.size());

        }

      }

      // Clear the map for further use.
      attribute_map.uniprotkb_id.clear();
      attribute_map.uniprotac_id.clear();
      attribute_map.attribute_map.clear();

      // Add the uniprot keys.
      attribute_map.uniprotkb_id = row_vector[VALUE_OFFSET_];
      attribute_map.uniprotac_id = row_vector[UNIPROTAC_OFFSET_];
      parser_state = PARSER_STATE::PROCESS_FIELD;
      // next record.
      continue;

    }

    if (parser_state == PARSER_STATE::INITIAL) {

      continue;

    }

    attribute_map.attribute_map.emplace(row_vector[FIELD_OFFSET_], row_vector[VALUE_OFFSET_]);

  }

  return true;

}


