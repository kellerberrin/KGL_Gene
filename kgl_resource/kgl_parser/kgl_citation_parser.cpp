//
// Created by kellerberrin on 26/7/21.
//

#include "kgl_citation_parser.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;


kgl::PMIDAlleleMap kgl::CitationResource::citationIndexedAlleles() const {

  PMIDAlleleMap pmid_indexed_alleles;

  for (auto const& [allele, pmid_vector] : citation_map_) {

    for (auto const& pmid : pmid_vector) {

      auto result = pmid_indexed_alleles.find(pmid);
      if (result == pmid_indexed_alleles.end()) {

        pmid_indexed_alleles.emplace(pmid, std::set<std::string>{allele});

      } else {

        auto& [pmid_key, allele_set] = *result;
        allele_set.insert(allele);

      }

    }

  }

  return pmid_indexed_alleles;

}


kgl::PMIDAlleleMap kgl::CitationResource::filteredCitationIndex(const std::set<std::string>& pmid_filter_set) const {

  PMIDAlleleMap pmid_indexed_alleles;

  for (auto const& [allele, pmid_vector] : citation_map_) {

    for (auto const& pmid : pmid_vector) {

      if (pmid_filter_set.contains(pmid)) {

        auto result = pmid_indexed_alleles.find(pmid);
        if (result == pmid_indexed_alleles.end()) {

          pmid_indexed_alleles.emplace(pmid, std::set<std::string>{allele});

        } else {

          auto& [pmid_key, allele_set] = *result;
          allele_set.insert(allele);

        }

      } // if in filter

    } // for pmid

  } // for allele

  return pmid_indexed_alleles;

}


kgl::DBCitationMap kgl::CitationResource::filteredAlleleIndexed(const std::set<std::string>& pmid_filter_set) const {

  DBCitationMap filtered_allele_map;

  auto filtered_pmid_map = filteredCitationIndex(pmid_filter_set);

  for (auto const& [pmid, allele_vector] : filtered_pmid_map) {

    for (auto const& allele : allele_vector) {

      auto result = filtered_allele_map.find(allele);
      if (result == filtered_allele_map.end()) {

        filtered_allele_map.emplace(allele, std::set<std::string>{pmid});

      } else {

        auto& [allele_key, pmid_set] = *result;
        pmid_set.insert(pmid);

      }

    }

  }

  return filtered_allele_map;

}



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

      citation_map_.emplace(rsid, std::set<std::string>{pmid_citation});

    } else {

      auto& [rsid, citation_set] = *find_result;
      citation_set.insert(pmid_citation);

    }

  }

  ExecEnv::log().info("Parsed Allele Index Citations; Alleles {} with Citations: {} from File: {}",
                      citation_map_.size() , parsed_record_ptr->getRowVector().size(), file_name);
  return true;

}
