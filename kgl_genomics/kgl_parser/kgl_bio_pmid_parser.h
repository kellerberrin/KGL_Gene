//
// Created by kellerberrin on 3/8/21.
//

#ifndef KGL_BIO_PMID_PARSER_H
#define KGL_BIO_PMID_PARSER_H


#include "kgl_resource_db.h"
#include "kgl_square_parser.h"



namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Gene nomenclature resource object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct BioPMIDRecord {

  std::string pmid_id;
  std::string bio_type;
  std::string bio_id;
  std::string bio_text;

};

using BioPMIDMap = std::multimap<std::string, BioPMIDRecord>;

class BioPMIDResource : public ResourceBase {

public:

  explicit BioPMIDResource(std::string identifier, BioPMIDMap&& disease_pmid_map, BioPMIDMap&& entrez_gene_map)
  : ResourceBase(std::move(identifier)), disease_pmid_map_(disease_pmid_map), entrez_pmid_map_(entrez_gene_map) {}
  ~BioPMIDResource() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::BIO_PMID; }

  [[nodiscard]] std::vector<std::string> entrezPMID(const std::string& entrez_id) const;   // In the form of an entrez gene identifier.
  [[nodiscard]] std::vector<std::string> diseaseMeSHPMID(const std::string& disease_mesh_id) const; // In the form of a mesh identifier, e.g. "MESH:D004487" (oedema).

  [[nodiscard]] const BioPMIDMap& entrezMap() const { return entrez_pmid_map_; }
  [[nodiscard]] const BioPMIDMap& diseaseMeSHMap() const { return disease_pmid_map_; }

private:

  const BioPMIDMap disease_pmid_map_;
  const BioPMIDMap entrez_pmid_map_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse a gene ident file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseBioPMID : public SquareTextParser {

public:

  ParseBioPMID() = default;
  ~ParseBioPMID() = default;

  [[nodiscard]] bool parseBioPMIDFile(const std::string& file_name);
  // Warning - for efficiency reasons the vectors are std::moved
  // The member vectors will be empty after the std::move.
  [[nodiscard]] BioPMIDMap&& moveEntrezPMIDMap() { return std::move(entrez_pmid_map_); }
  [[nodiscard]] BioPMIDMap&& moveDiseasePMIDMap() { return std::move(disease_pmid_map_); }

private:

  BioPMIDMap  entrez_pmid_map_;
  BioPMIDMap  disease_pmid_map_;

  constexpr static const size_t MINIMUM_ROW_COUNT_{1};
  constexpr static const size_t COLUMN_COUNT_{5};

  constexpr static const char* DISEASE_TAG{"Disease"};
  constexpr static const char* ENTREZ_GENE_TAG{"Gene"};
  // Field record offsets.
  constexpr static const size_t PMID_OFFSET_{0};
  constexpr static const size_t BIO_TYPE_OFFSET_{1};
  constexpr static const size_t BIO_ID_OFFSET_{2};
  constexpr static const size_t BIO_TEXT_OFFSET_{3};
  constexpr static const size_t BIO_RESOURCE_OFFSET_{4};

};


} // Namespace


#endif //KGL_BIO_PMID_PARSER_H
