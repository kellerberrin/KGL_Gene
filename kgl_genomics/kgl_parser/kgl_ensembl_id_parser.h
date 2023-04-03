//
// Created by kellerberrin on 16/3/21.
//

#ifndef KGL_ENSEMBL_ID_PARSER_H
#define KGL_ENSEMBL_ID_PARSER_H

#include "kgl_properties_resource.h"
#include "kgl_square_parser.h"


namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Gene nomenclature resource object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct GeneIDSynonyms {

  std::string HGNC_id;
  std::string Canonical;
  std::string ensembl_id;

};
using GeneSynonymVector = std::vector<GeneIDSynonyms>;
using HGNCEnsemblMap = std::map<std::string, std::string>;

class EnsemblHGNCResource : public ResourceBase {

public:

  explicit EnsemblHGNCResource(std::string identifier, GeneSynonymVector synonym_vector)
    : ResourceBase(ResourceProperties::GENE_NOMENCLATURE_RESOURCE_ID_, std::move(identifier)), synonym_vector_(std::move(synonym_vector)) {

    IndexHGNC();

  }
  ~EnsemblHGNCResource() override = default;

  [[nodiscard]] std::string HGNCToEnsembl(const std::string& hgnc_id) const;

private:

  const GeneSynonymVector synonym_vector_;
  HGNCEnsemblMap hgnc_emsembl_map_;
  void IndexHGNC();

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse a gene ident file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseGeneIdents : public SquareTextParser {

public:

  ParseGeneIdents() = default;
  ~ParseGeneIdents() = default;

  [[nodiscard]] bool parseIdentFile(const std::string& file_name);
  [[nodiscard]] const GeneSynonymVector& getSynonymVector() const { return synonym_vector_; }

private:

  GeneSynonymVector  synonym_vector_;

  constexpr static const size_t MINIMUM_ROW_COUNT_{1};
  constexpr static const size_t COLUMN_COUNT_{3};

  constexpr static const size_t HGNC_OFFSET_{2};
  constexpr static const size_t CANONICAL_OFFSET_{1};
  constexpr static const char* CANONICAL_VALUE_{"1"};
  constexpr static const size_t ENSEMBL_OFFSET_{0};

};


} // namespace


#endif //KGL_ENSEMBL_ID_PARSER_H
