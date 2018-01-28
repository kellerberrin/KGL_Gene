//
// Created by kellerberrin on 26/01/18.
//

#include <fstream>

#include "kgl_exec_env.h"
#include "kgl_gaf_parser.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


namespace bt = boost;
namespace kgl = kellerberrin::genome;



bool kgl::GeneOntology::readGafFile(const std::string &file_name) {


  ExecEnv::log().info("Reading Gaf file: {}", file_name);

  std::ifstream gaf_file;

  // Open input file.

  gaf_file.open(file_name);

  if (not gaf_file.good()) {

    ExecEnv::log().critical("I/O error; could not open Gaf file: {}", file_name);

  }

  try {

    long counter = 0;

    while (true) {

      std::string record_str;

      if (std::getline(gaf_file, record_str).eof()) break;

      if ((record_str)[0] == '!') continue;   // ignore header records.

      parseGafRecord(record_str);

      ++counter;

      if (counter % REPORT_INCREMENT_ == 0) {

        ExecEnv::log().info("Read: {} Gaf records", counter);

      }

    }

    gaf_file.close();

    ExecEnv::log().info("Processed: {} Gaf records", counter);

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("Gaf file: {}, unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  return true;

}


void kgl::GeneOntology::parseGafRecord(const std::string& record_str) {

  std::vector<std::string> field_vec;
  bt::char_separator<char> item_key_sep("\t");
  bt::tokenizer<bt::char_separator<char>> tokenize_item(record_str, item_key_sep);
  for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

    field_vec.push_back(*iter_item);

  }

}

