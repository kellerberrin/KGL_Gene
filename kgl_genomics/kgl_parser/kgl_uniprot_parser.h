//
// Created by kellerberrin on 2/7/21.
//

#ifndef KGL_UNIPROT_PARSER_H
#define KGL_UNIPROT_PARSER_H



#include "kgl_resource_db.h"
#include "kgl_square_parser.h"


namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Uniprot gene info resource.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// key= field_id, value = field_value
struct UniprotAttributeMap {

  std::string uniprotkb_id;     // Uniprot knowledge base id, unique index
  std::string uniprotac_id;     // Uniprot accession id, non-unique index;
  std::multimap<std::string, std::string> attribute_map; // {field_id, field_value} pair (non unique).

};
// key = uniprot_id , value = field_map (above).
using UniprotIDMap = std::map<std::string, UniprotAttributeMap>;

class UniprotResource : public ResourceBase {

public:

  explicit UniprotResource(std::string identifier, UniprotIDMap&& uniprot_map)
  : ResourceBase(std::move(identifier)), uniprot_map_(uniprot_map) {

    createIndexes();

  }
  ~UniprotResource() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::GENE_NOMENCLATURE; }

  [[nodiscard]] std::vector<std::string> symbolToUniprot(const std::string& symbol) const;
  [[nodiscard]] std::vector<std::string> HGNCToEnsembl(const std::string& hgnc) const { return lookupInfo( hgnc, hgnc_index_, ENSEMBL_FIELD); }
  [[nodiscard]] std::vector<std::string> HGNCToUniprot(const std::string& hgnc) const;
  [[nodiscard]] std::vector<std::string> uniprotToEnsembl(const std::string& hgnc) const;

  const static constexpr char* UNIPROTKB_ID{"UniProtKB-ID"};
  const static constexpr char* ENTREZ_GENE{"GeneID"};
  const static constexpr char* GENE_NAME{"Gene_Name"};
  const static constexpr char* GENE_SYNONYM{"Gene_Synonym"};
  const static constexpr char* ENSEMBL_FIELD{"Ensembl"};
  const static constexpr char* HGNC_FIELD{"HGNC"};

private:

  const UniprotIDMap uniprot_map_;
  std::multimap<std::string, std::string> uniprot_index_;
  std::multimap<std::string, std::string> entrez_index_;
  std::multimap<std::string, std::string> symbol_index_;
  std::multimap<std::string, std::string> synonym_index_;
  std::multimap<std::string, std::string> ensembl_index_;
  std::multimap<std::string, std::string> hgnc_index_;

  void createIndexes();

  [[nodiscard]] std::vector<std::string> lookupInfo(const std::string& lookup_value,
                                                    const std::multimap<std::string, std::string>& lookup_map,
                                                    const std::string& lookup_type) const;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse the Uniprot Info Resource file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseUniprotId : public SquareTextParser {

public:

  ParseUniprotId() = default;
  ~ParseUniprotId() = default;

  [[nodiscard]] bool parseUniprotFile(const std::string& file_name);
  [[nodiscard]] UniprotIDMap&& getUniproResource() { return std::move(uniprot_info_); };

private:

  UniprotIDMap  uniprot_info_;

  constexpr static const size_t COLUMN_COUNT_{3};
  constexpr static const size_t UNIPROTAC_OFFSET_{0};
  constexpr static const size_t FIELD_OFFSET_{1};
  constexpr static const size_t VALUE_OFFSET_{2};

};


} // namespace


#endif //KGL_UNIPROT_PARSER_H
