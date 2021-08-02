//
// Created by kellerberrin on 26/7/21.
//

#ifndef KGL_CITATION_PARSER_H
#define KGL_CITATION_PARSER_H


#include "kgl_resource_db.h"
#include "kgl_square_parser.h"
#include "kgl_json_parser.h"


namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Allele  Citation resource object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CitationResource : public ResourceBase {

public:

  CitationResource(std::string identifier, DBCitationMap citation_map) : ResourceBase(std::move(identifier)), citation_map_(std::move(citation_map))  {}
  ~CitationResource() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::ALLELE_CITATION; }

  [[nodiscard]] const DBCitationMap& citationMap() const { return citation_map_; }

private:

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
