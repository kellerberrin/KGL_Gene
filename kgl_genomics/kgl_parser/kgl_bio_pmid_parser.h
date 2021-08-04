//
// Created by kellerberrin on 3/8/21.
//

#ifndef KGL_BIO_PMID_PARSER_H
#define KGL_BIO_PMID_PARSER_H


#include "kgl_resource_db.h"
#include "kgl_data_file_type.h"

#include <set>


namespace kellerberrin::genome {   //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// PMID Pubmed publication resource object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct BioPMIDRecord {

  BioPMIDRecord() = default;
  BioPMIDRecord(const BioPMIDRecord& record);
  BioPMIDRecord(  std::string&& pmid_id,
                  std::string&& bio_type,
                  std::string&& bio_id,
                  std::string&& bio_text) noexcept;
  BioPMIDRecord(BioPMIDRecord&& record) noexcept;
  BioPMIDRecord& operator=(BioPMIDRecord&& record) noexcept;

  std::string pmid_id;
  std::string bio_type;
  std::string bio_id;
  std::string bio_text;


};

using BioPMIDMap = std::multimap<std::string, BioPMIDRecord>;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Holds the PMID data indexed by Gene Entrez id (all species) and Disease Mesh id (all diseases).
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BioPMIDMaps {

public:

  BioPMIDMaps(BioPMIDMap&& disease_pmid_map, BioPMIDMap&& entrez_gene_map) : disease_pmid_map_(disease_pmid_map), entrez_pmid_map_(entrez_gene_map) {}
  ~BioPMIDMaps() = default;

  // Only want unique pmid identifiers.
  [[nodiscard]] std::set<std::string> entrezPMID(const std::string& entrez_id) const;   // In the form of an entrez gene identifier.
  [[nodiscard]] std::set<std::string> diseaseMeSHPMID(const std::string& disease_mesh_id) const; // In the form of a mesh identifier, e.g. "MESH:D004487" (oedema).

  [[nodiscard]] const BioPMIDMap& entrezMap() const { return entrez_pmid_map_; }
  [[nodiscard]] const BioPMIDMap& diseaseMeSHMap() const { return disease_pmid_map_; }

private:

  BioPMIDMap disease_pmid_map_;
  BioPMIDMap entrez_pmid_map_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PMID Resource object. Generally the PMID data is too large to be loaded as a resource (~100GB).
// This object makes resource functionality available if available memory permits.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BioPMIDResource : public ResourceBase {

public:

  explicit BioPMIDResource(std::string identifier, BioPMIDMap&& disease_pmid_map, BioPMIDMap&& entrez_gene_map)
  : ResourceBase(std::move(identifier)), pmid_maps_(std::move(disease_pmid_map), std::move(entrez_gene_map)) {}
  ~BioPMIDResource() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::BIO_PMID; }

  // Only want unique pmid identifiers.
  [[nodiscard]] std::set<std::string> entrezPMID(const std::string& entrez_id) const { return pmid_maps_.entrezPMID(entrez_id); }
  [[nodiscard]] std::set<std::string> diseaseMeSHPMID(const std::string& disease_mesh_id) const { return pmid_maps_.diseaseMeSHPMID(disease_mesh_id); }

  [[nodiscard]] const BioPMIDMap& entrezMap() const { return pmid_maps_.entrezMap(); }
  [[nodiscard]] const BioPMIDMap& diseaseMeSHMap() const { return pmid_maps_.diseaseMeSHMap(); }

private:

  BioPMIDMaps pmid_maps_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PMID File object. Loaded into the requesting package as a file object.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BioPMIDFileData : public DataDB {

public:

  BioPMIDFileData(const std::string& file_name, BioPMIDMap&& disease_pmid_map, BioPMIDMap&& entrez_gene_map)
    : DataDB(DataSourceEnum::BioPMID), file_name_(file_name), pmid_maps_(std::move(disease_pmid_map), std::move(entrez_gene_map))  {}
  ~BioPMIDFileData() override = default;

  [[nodiscard]] const std::string& fileId() const override { return file_name_; }

  // Only want unique pmid identifiers.
  [[nodiscard]] std::set<std::string> entrezPMID(const std::string& entrez_id) const { return pmid_maps_.entrezPMID(entrez_id); }
  [[nodiscard]] std::set<std::string> diseaseMeSHPMID(const std::string& disease_mesh_id) const { return pmid_maps_.diseaseMeSHPMID(disease_mesh_id); }

  [[nodiscard]] const BioPMIDMap& entrezMap() const { return pmid_maps_.entrezMap(); }
  [[nodiscard]] const BioPMIDMap& diseaseMeSHMap() const { return pmid_maps_.diseaseMeSHMap(); }

private:

  std::string file_name_;
  BioPMIDMaps pmid_maps_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse a pmid bio resource.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseBioPMID {

public:

  ParseBioPMID() = default;
  ~ParseBioPMID() = default;

  [[nodiscard]] bool parseBioPMIDRecords(const std::string& file_name);
  // Warning - for efficiency reasons the vectors are std::moved
  // The member vectors will be empty after the std::move.
  [[nodiscard]] BioPMIDMap&& moveEntrezPMIDMap() { return std::move(entrez_pmid_map_); }
  [[nodiscard]] BioPMIDMap&& moveDiseasePMIDMap() { return std::move(disease_pmid_map_); }

private:

  BioPMIDMap  entrez_pmid_map_;
  BioPMIDMap  disease_pmid_map_;

  constexpr static const size_t COLUMN_COUNT_{5};
  constexpr static const size_t REPORT_INTERVAL_{1000000};

  constexpr static const char COMMENT_{'#'}; // If first character then line ignored.
  constexpr static const char DELIMITER_{'\t'}; // Delimits the fields on each line.

  constexpr static const char* DISEASE_TAG{"Disease"};
  constexpr static const char* ENTREZ_GENE_TAG{"Gene"};
  // Field record offsets.
  constexpr static const size_t PMID_OFFSET_{0};
  constexpr static const size_t BIO_TYPE_OFFSET_{1};
  constexpr static const size_t BIO_ID_OFFSET_{2};
  constexpr static const size_t BIO_TEXT_OFFSET_{3};
  constexpr static const size_t BIO_RESOURCE_OFFSET_{4};

  bool parseFields(const std::vector<std::string_view>& field_views);

};


} // Namespace


#endif //KGL_BIO_PMID_PARSER_H
