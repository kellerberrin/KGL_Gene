//
// Created by kellerberrin on 23/7/21.
//

#ifndef KGL_VARIANT_INFO_JSON_H
#define KGL_VARIANT_INFO_JSON_H


#include "kgl_genome_types.h"
#include "kgl_data_file_impl.h"
#include "kgl_variant_db_type.h"

#include "kel_bound_queue.h"

#include <memory>
#include <string>
#include <vector>
#include <map>


namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This DB object is passed back to the requesting package for further processing (typically written to a file).
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using DBCitationMap = std::map<GenomeId_t, std::vector<std::string>>;

class DBCitation : public DataDB {

public:

  DBCitation(DataSourceEnum data_source, const std::string& file_name) : DataDB(data_source), file_name_(file_name) {}
  ~DBCitation() override = default;

  [[nodiscard]] const DBCitationMap& citationMap() const { return citation_map_; }
  [[nodiscard]] bool insertCitations(const std::string& allele_rsid, const std::vector<std::string>& pmid_citations);
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

  // Begin reading IO records, spawns threads.
  [[nodiscard]] bool commenceJSONIO(const std::string& json_file_name);
  // Parse the Json file for PMID citations.
  void parseJson(const std::shared_ptr<DBCitation>& db_citation_ptr);
  [[nodiscard]] const std::string& getFileName() const { return file_data_.fileName(); }

private:

  // The upstream queue of line records.
  FileDataIO file_data_;

  static const constexpr size_t REPORT_INTERVAL{100000}; // Parser progress messages.

};



} // end namespace



#endif // KGL_VARIANT_INFO_JSON_H
