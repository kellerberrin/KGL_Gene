//
// Created by kellerberrin on 23/7/21.
//

#ifndef KGL_VARIANT_INFO_JSON_H
#define KGL_VARIANT_INFO_JSON_H


#include "kgl_genome_types.h"
#include "kgl_data_file_impl.h"
#include "kgl_data_file_type.h"

#include "../../kel_thread/kel_queue_tidal.h"

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <set>


namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This DB object is passed back to the requesting package for further processing.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The map is sorted by variant key=rsid, value=set of pmid (Pubmed article identifiers)
using DBCitationMap = std::map<std::string, std::set<std::string>>;

class DBCitation : public DataDB {

public:

  DBCitation(DataSourceEnum data_source, const std::string& file_name) : DataDB(data_source), file_name_(file_name) {}
  ~DBCitation() override = default;

  [[nodiscard]] const DBCitationMap& citationMap() const { return citation_map_; }
  [[nodiscard]] DBCitationMap& citationMap() { return citation_map_; }
  [[nodiscard]] const std::string& fileId() const override { return file_name_; }

private:

  const std::string file_name_;
  DBCitationMap citation_map_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parses a Json file for citations, currently uses the 3rd party library 'RapidJson'. Any Json parser should be OK.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class JSONInfoParser {

public:

  JSONInfoParser() = default;
  ~JSONInfoParser() = default;

  [[nodiscard]] bool parseFile(const std::string& json_file_name, const std::shared_ptr<DBCitation>& db_citation_ptr);
  [[nodiscard]] bool parseFile(const std::string& json_file_name, DBCitationMap& citation_map);

  [[nodiscard]] const std::string& getFileName() const { return file_data_.fileName(); }

private:

  // The upstream queue of line records.
  FileDataIO file_data_;

  static const constexpr size_t REPORT_INTERVAL_{100000}; // Parser progress messages.

  // Begin reading IO records, spawns threads.
  [[nodiscard]] bool commenceJSONIO(const std::string& json_file_name);
  // Parse the Json file for PMID citations.
  void parseJson(DBCitationMap& citation_map);

};



} // end namespace



#endif // KGL_VARIANT_INFO_JSON_H
