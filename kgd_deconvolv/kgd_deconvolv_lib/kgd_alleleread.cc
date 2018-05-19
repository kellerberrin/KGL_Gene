//
// Created by kellerberrin on 19/05/18.
//

#include "kgd_alleleread.h"
#include "kgd_utility.h"

#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

namespace bt = boost;
namespace kgd = kellerberrin::deconvolv;


bool kgd::AlleleReader::parseFile(const std::string& file_name, const std::string& delimiters) {

  std::ifstream data_file;
  size_t counter = 0;

  // Open input file.

  data_file.open(file_name);

  if (not data_file.good()) {

    ExecEnv::log().critical("alleleReader; I/O error; could not open data file: {}", file_name);

  }

  try {

    std::string record_str;
    std::vector<std::string> text_fields;
    bt::char_separator<char> item_key_sep(delimiters.c_str());

    while (not std::getline(data_file, record_str).eof()) {

      if(record_str.size() == 0) continue;

      if (record_str[0] == COMMENT_CHARACTER_) continue;   // ignore header records.

      bt::tokenizer<bt::char_separator<char>> tokenize_item(record_str, item_key_sep);

      for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

        text_fields.push_back(Utility::trimWhiteSpace(*iter_item));

      }

      if (not parseDataLine(text_fields)) {

        ExecEnv::log().error("AlleleReader, Problem parsing file: {}, line record: {}", file_name, record_str);

      }

      ++counter;

    }

    data_file.close();

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("AlleleReader, File: {}, unexpected I/O exception: {}", file_name, e.what());

  }

  ExecEnv::log().info("AlelleReader, Parsed: {} lines from aux file: {}", counter, file_name);

  return true;

}


bool kgd::AlleleReader::parseDataLine(const std::vector<std::string>& text_fields) {


  return true;

}

