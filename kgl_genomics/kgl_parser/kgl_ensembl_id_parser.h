//
// Created by kellerberrin on 16/3/21.
//

#ifndef KGL_ENSEMBL_ID_PARSER_H
#define KGL_ENSEMBL_ID_PARSER_H

#include "kgl_square_parser.h"


namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse a gene ident file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct GeneIDSynonyms {

  std::string HGNC_id;
  std::string ensembl_id;

};
using GeneSynonymVector = std::vector<GeneIDSynonyms>;


class ParseGeneIdents : public SquareTextParser {

public:

  ParseGeneIdents() = default;
  ~ParseGeneIdents() = default;

  [[nodiscard]] bool parseIdentFile(const std::string& file_name);
  [[nodiscard]] const GeneSynonymVector& getSynonymVector() const { return synonym_vector_; }

private:

  GeneSynonymVector  synonym_vector_;

  constexpr static const size_t MINIMUM_ROW_COUNT_{1};
  constexpr static const size_t COLUMN_COUNT_{2};

  constexpr static const size_t HGNC_OFFSET_{1};
  constexpr static const size_t ENSEMBL_OFFSET_{0};

};


} // namespace


#endif //KGL_ENSEMBL_ID_PARSER_H
