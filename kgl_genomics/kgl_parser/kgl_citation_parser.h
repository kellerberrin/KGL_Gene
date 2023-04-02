//
// Created by kellerberrin on 26/7/21.
//

#ifndef KGL_CITATION_PARSER_H
#define KGL_CITATION_PARSER_H


#include "kgl_runtime_resource.h"
#include "kgl_square_parser.h"
#include "kgl_json_parser.h"

#include <map>
#include <set>

namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Allele  Citation resource object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Key is the Pubmed publication pmid, value is a set of unique alleles as 'rsXXX' identifiers.
using PMIDAlleleMap = std::map<std::string, std::set<std::string>>;

class CitationResource : public ResourceBase {

public:

  CitationResource(std::string identifier, DBCitationMap citation_map) : ResourceBase(std::move(identifier)), citation_map_(std::move(citation_map))  {}
  ~CitationResource() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::ALLELE_CITATION; }

  // Indexed by allele 'rsXXX' code, value is a vector of Pubmed publication pmids.
  [[nodiscard]] const DBCitationMap& alleleIndexedCitations() const { return citation_map_; }
  // Filtered by allele 'rsXXX' code, value is a vector of Pubmed publication pmids present in the filtered set.
  [[nodiscard]] DBCitationMap filteredAlleleIndexed(const std::set<std::string>& pmid_filter_set) const;
  // Citation indexed alleles. Key is pmid, value is a vector of alleles 'rsXXX' code.
  [[nodiscard]] PMIDAlleleMap citationIndexedAlleles() const;
  // Filtered citation indexed alleles, only pmids/alleles in the filter set are returned.
  [[nodiscard]] PMIDAlleleMap filteredCitationIndex(const std::set<std::string>& pmid_filter_set) const;


private:

  // Indexed by allele 'rsXXX' code, value is a unique set of Pubmed publication pmids.
  const DBCitationMap citation_map_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse a gene ident file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseCitations : public SquareTextParser {

public:

  ParseCitations() = default;
  ~ParseCitations() = default;

  [[nodiscard]] bool parseCitationFile(const std::string& file_name);
  [[nodiscard]] const DBCitationMap& getCitationMap() const { return citation_map_; }

private:

  DBCitationMap citation_map_;

  constexpr static const size_t MINIMUM_ROW_COUNT_{1};
  constexpr static const size_t COLUMN_COUNT_{2};

  constexpr static const size_t RSID_OFFSET_{0};
  constexpr static const size_t PMID_OFFSET_{1};

};


} // namespace


#endif // KGL_CITATION_PARSER_H
