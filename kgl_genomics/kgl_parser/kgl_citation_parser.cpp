//
// Created by kellerberrin on 26/7/21.
//

#include "kgl_citation_parser.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;


[[nodiscard]] bool kgl::ParseCitations::parseCitationFile(const std::string &file_name) {


  auto parsed_record_ptr = parseFlatFile(file_name, SquareTextParser::DELIMITER_CSV);

  if (parsed_record_ptr->getRowVector().size() < MINIMUM_ROW_COUNT_) {

    ExecEnv::log().error("ParseCitations::parseCitationFile; Row count: {} for file: {} is below minimum",
                         parsed_record_ptr->getRowVector().size(), file_name);
    return false;

  }

  if (not parsed_record_ptr->checkRowSize(COLUMN_COUNT_)) {

    ExecEnv::log().error("ParseCitations::parseCitationFile; Not all rows have expected column count: {} for file: {}",
                         COLUMN_COUNT_, file_name);
    return false;

  }

  // No header with this file.
  for (const auto& row_vector :  parsed_record_ptr->getRowVector()) {

    auto rsid = Utility::trimEndWhiteSpace(row_vector[RSID_OFFSET_]);
    auto pmid_citation = Utility::trimEndWhiteSpace(row_vector[PMID_OFFSET_]);

    auto find_result = citation_map_.find(rsid);
    if (find_result == citation_map_.end()) {

      citation_map_.emplace(rsid, std::vector<std::string>{pmid_citation});

    } else {

      auto& [rsid, citation_vector] = *find_result;
      citation_vector.push_back(pmid_citation);

    }

  }

  ExecEnv::log().info("Parsed Allele Index Citations; Alleles {} with Citations: {} from File: {}",
                      citation_map_.size() , parsed_record_ptr->getRowVector().size(), file_name);
  return true;

}
