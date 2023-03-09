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

  file_name_ = json_file_name;
  if (not file_data_.open(json_file_name)) {

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

    IOLineRecord line_record = file_data_.readLine();
    if (line_record.EOFRecord()) { // check for EOF condition.

      break;

    }

    auto const [line_count, json_line] = line_record.getLineData();

    chars_processed += json_line.size();
    ++lines_processed;

    if (lines_processed % REPORT_INTERVAL_ == 0) {

      ExecEnv::log().info("JSONInfoParser processed line count: {}, text bytes: {}, file: {}",
                          line_count, chars_processed, getFileName());

    }

    rapidjson::Document document;
    document.Parse(json_line.data(), json_line.size());

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
        std::set<std::string> cite_string_set;
        for (size_t idx = 0; idx < cite_array.Size(); ++idx) {

          rapidjson::Value& pmid = cite_array[idx];
          size_t pmid_cite = pmid.GetInt64();
          cite_string_set.insert(std::to_string(pmid_cite));

        }

        auto [iter, result] = citation_map.try_emplace(rsid, cite_string_set);
        if (not result) {

          ExecEnv::log().error("JSONInfoParser::parseJson; problem processing citations for (duplicate) rsid: {}", rsid);

        }

      }

    }

  } // while

  ExecEnv::log().info("JSONInfoParser::dequeueJSONline; records/lines: {}, chars: {}", lines_processed, chars_processed);
  ExecEnv::log().info("JSONInfoParser::dequeueJSONline; cited alleles: {}, total citations: {}", cited_alleles, total_citations);

}
