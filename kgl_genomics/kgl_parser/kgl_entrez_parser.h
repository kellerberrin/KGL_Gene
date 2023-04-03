//
// Created by kellerberrin on 2/8/21.
//

#ifndef KGL_ENTREZ_PARSER_H
#define KGL_ENTREZ_PARSER_H

#include "kgl_properties_resource.h"
#include "kgl_square_parser.h"


namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Gene nomenclature resource object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct EntrezRecord {

  std::string entrez_id;
  std::string symbol_id;
  std::string description;

};
using EntrezVector = std::vector<EntrezRecord>;
using EntrezMap = std::map<std::string, EntrezRecord>;

class EntrezResource : public ResourceBase {

public:

  explicit EntrezResource(std::string identifier, EntrezVector entrez_vector)
    : ResourceBase(ResourceProperties::ENTREZ_RESOURCE_ID_, std::move(identifier)),
  entrez_vector_(std::move(entrez_vector)) {

    IndexSymbol();
    IndexEntrez();

  }
  ~EntrezResource() override = default;

  [[nodiscard]] std::string symbolToEntrez(const std::string& symbol_id) const;
  [[nodiscard]] std::string entrezToSymbol(const std::string& entrez_id) const;

private:

  const EntrezVector entrez_vector_;
  EntrezMap symbol_map_;
  EntrezMap entrez_map_;

  void IndexSymbol();
  void IndexEntrez();

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse a gene ident file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseEntrez : public SquareTextParser {

public:

  ParseEntrez() = default;
  ~ParseEntrez() = default;

  [[nodiscard]] bool parseEntrezFile(const std::string& file_name);
  [[nodiscard]] const EntrezVector& getEntrezVector() const { return entrez_vector_; }

private:

  EntrezVector  entrez_vector_;

  constexpr static const size_t MINIMUM_ROW_COUNT_{1};
  constexpr static const size_t COLUMN_COUNT_{16};

  // Field record offsets.
  constexpr static const size_t TAXON_ID_OFFSET_{0};
  constexpr static const size_t ENTREZ_OFFSET_{1};
  constexpr static const size_t SYMBOL_OFFSET_{2};
  constexpr static const size_t LOCUSTAG_OFFSET_{3};
  constexpr static const size_t SYNONYMS_OFFSET_{4};
  constexpr static const size_t CHROMOSOME_OFFSET_{5};
  constexpr static const size_t DBXREFS_OFFSET_{6};
  constexpr static const size_t MAP_LOCATION_OFFSET_{7};
  constexpr static const size_t DESCRIPTION_OFFSET_{8};
  constexpr static const size_t TYPE_OF_GENE_OFFSET_{9};
  constexpr static const size_t SYMBOL_FROM_NOMENCLATURE_AUTHORITY_OFFSET_{10};
  constexpr static const size_t FULL_NAME_FROM_NOMENCLATURE_AUTHORITY_OFFSET_{11};
  constexpr static const size_t NOMENCLATURE_STATUS_OFFSET_{12};
  constexpr static const size_t OTHER_DESIGNATIONS_OFFSET_{13};
  constexpr static const size_t MODIFICATION_DATE_OFFSET_{14};
  constexpr static const size_t FEATURE_TYPE_OFFSET_{15};

};


} // namespace


#endif //KGL_ENTREZ_PARSER_H
