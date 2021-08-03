//
// Created by kellerberrin on 23/7/21.
//

#include "kgl_json_parser.h"

#include "rapidjson/document.h"
#include "rapidjson/error/en.h"
#include "rapidjson/writer.h"

#include <sstream>
#include <memory>


namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::JSONInfoParser::parseFile(const std::string& json_file_name, const std::shared_ptr<DBCitation>& db_citation_ptr) {

  return parseFile(json_file_name, db_citation_ptr->citationMap());

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::JSONInfoParser::parseFile(const std::string& json_file_name, DBCitationMap& citation_map) {

  if (commenceJSONIO(json_file_name)) {

    parseJson(citation_map);
    ExecEnv::log().info("JSON File: {} has alleles with PMID citations: {}", json_file_name, citation_map.size());

  } else {

    return false;

  }

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::JSONInfoParser::commenceJSONIO(const std::string& json_file_name) {

  if (not file_data_.commenceIO(json_file_name)) {

    ExecEnv::log().error("JSONInfoParser::commenceJSONIO; problem processing JSON file: {}", json_file_name);
    return false;

  }

  return true;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dequeue a line record and incrementally call the JSON parser

void kgl::JSONInfoParser::parseJson(DBCitationMap& citation_map) {

  size_t chars_processed{0};
  size_t lines_processed{0};
  size_t cited_alleles{0};
  size_t total_citations{0};
  std::string concatonated;
  std::vector<std::unique_ptr<std::string>> record_lines;
  const std::string citation_key{"citations"};
  const std::string ref_snp_key{"refsnp_id"};

  while (true) {

    IOLineRecord line_record = file_data_.readIORecord();
    if (not line_record) { // check for EOF condition.

      // Push the EOF marker back on the queue for other threads.
      file_data_.enqueueEOF();
      break;

    }

    auto& [line_count, line_string_ptr] = line_record.value();

    size_t line_length = line_string_ptr->size();
    chars_processed += line_length;
    ++lines_processed;

    if (lines_processed % REPORT_INTERVAL == 0) {

      ExecEnv::log().info("JSONInfoParser processed line count: {}, text bytes: {}, file: {}",
                          line_count, chars_processed, getFileName());

    }

    rapidjson::Document document;
    document.Parse(line_string_ptr->data(), line_string_ptr->size());

    if (document.HasParseError()) {

      ExecEnv::log().warn("JSONInfoParser::dequeueJSONline; error parsing Json line at offset: {}, error text: {}",
                          document.GetErrorOffset(), GetParseError_En(document.GetParseError()));
      continue;

    }

    rapidjson::Value& citations = document[citation_key.c_str()];
    if (citations.IsArray()) {

      auto cite_array = citations.GetArray();
      if (cite_array.Size() > 0) {

        ++cited_alleles;
        total_citations += cite_array.Size();
        rapidjson::Value& ref_snp = document[ref_snp_key.c_str()];

        std::string rsid{"rs"};
        rsid += ref_snp.GetString();
        std::vector<std::string> cite_string_array;
        for (size_t idx = 0; idx < cite_array.Size(); ++idx) {

          rapidjson::Value& pmid = cite_array[idx];
          size_t pmid_cite = pmid.GetInt64();
          cite_string_array.emplace_back(std::to_string(pmid_cite));

        }

        auto [iter, result] = citation_map.try_emplace(rsid, cite_string_array);
        if (not result) {

          ExecEnv::log().error("JSONInfoParser::parseJson; problem processing citations for (duplicate) rsid: {}", rsid);

        }

      }

    }

  } // while

  ExecEnv::log().info("JSONInfoParser::dequeueJSONline; records/lines: {}, chars: {}", lines_processed, chars_processed);
  ExecEnv::log().info("JSONInfoParser::dequeueJSONline; cited alleles: {}, total citations: {}", cited_alleles, total_citations);

}
