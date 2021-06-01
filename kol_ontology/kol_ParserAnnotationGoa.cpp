//
// Created by kellerberrin on 26/5/21.
//
#include "kol_ParserAnnotationGoa.h"

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <iostream>
#include <boost/tokenizer.hpp>


namespace kol = kellerberrin::ontology;


std::shared_ptr<kol::AnnotationData> kol::ParserAnnotationGoa::parseAnnotationFile(const std::string &file_name) const {

  std::vector<std::shared_ptr<const GAFRecord>> gaf_records;
  // Check that the supplied policy is valid.
  if (not _policy.isValid()) {

    ExecEnv::log().error("ParserAnnotationGoa::parseAnnotationFile; invalid annotation policy for file: {}", file_name);
    return std::make_shared<AnnotationData>(gaf_records);

  }

  std::ifstream anno_file(file_name);
  if (not anno_file.good()) {

    ExecEnv::log().error("ParserAnnotationGoa::parseAnnotationFile; problem opening annotation file: {}", file_name);
    return  std::make_shared<AnnotationData>(gaf_records);

  }

  std::string line;
  size_t line_count{0};
  while (not std::getline(anno_file, line).eof()) {

    ++line_count;
    if (line[0] == COMMENT_CHAR) {

      continue;

    }

    auto view_vector = Utility::view_tokenizer(line, TAB_FIELD_DELIMITER);
    if (view_vector.size() != EXPECTED_FIELD_COUNT) {

      ExecEnv::log().error("ParserAnnotationGoa::parseAnnotationFile; line: {}, expected fields: {}, found: {}, file: {}",
                           line_count, EXPECTED_FIELD_COUNT, view_vector.size(), file_name);
      continue;

    }

    if (not view_vector[FIELD_OFFSET_QUALIFIER].empty()) {

      auto qualifier = Utility::toupper(std::string(view_vector[FIELD_OFFSET_QUALIFIER]));
      if (qualifier.find(NOT_QUALIFIER) != std::string::npos) {

        continue;

      }

    }

    std::shared_ptr<GAFRecord> record_ptr(std::make_shared<GAFRecord>());
    if (not record_ptr->parseGafRecord(line)) {

      ExecEnv::log().error("AnnotationData::addGAFRecord; error parsing Gaf record from line: {}", line);

    } else {

      gaf_records.push_back(record_ptr);

    }

  }

  auto filtered_records = AnnotationData::filterGAFRecords(_policy, gaf_records);
  std::shared_ptr<AnnotationData> annoData(std::make_shared<AnnotationData>(filtered_records));

  return annoData;

}


//! A method for checking if a file exists and is formatted correctly.
/*!
  This function checks that the file exists and its format can be recognized.
*/

bool kol::ParserAnnotationGoa::isFileGood(const std::string &fileName) const {

  std::ifstream in(fileName.c_str());

  if (!in.good()) {

    return false;

  }

  //Tokenizer type
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  //declares a tab separator variable
  boost::char_separator<char> tab_sep("\t", "", boost::keep_empty_tokens);
  tokenizer::iterator it;

  std::size_t count = 0;
  std::string line;

  while (in.good() && count < 5) {
    //get next line in 'in' file stream
    std::getline(in, line);

    //split line
    tokenizer tokens(line, tab_sep);
    //set iterator to first token
    it = tokens.begin();

    //always check if empty (must have in linux)
    if (it == tokens.end()) { continue; }

    //skip comments if line is not empty and starts with !
    if (line.at(0) == '!') { continue; }


    std::string database, geneName, qualifierStr, goString, evidenceCode, ontology;
    std::size_t i = 0;
    for (; it != tokens.end(); ++it) {
      switch (i) {
        case 0:
          database = *it;
          if (database.size() == 0) { return false; }
          break;
        case 1:
          geneName = *it;
          if (geneName.size() == 0) { return false; }
          break;
        case 3:
          qualifierStr = *it;
          break;
        case 4:
          goString = *it;
          if (goString.size() == 0) { return false; }           // disallow empty go
          if (goString.substr(0, 3) != "GO:") { return false; } // disallow bad go term
          break;
        case 6:
          evidenceCode = *it;
          if (evidenceCode.size() == 0) { return false; }
          if (GO::evidenceStringToCode(evidenceCode) == GO::EvidenceCode::ECODE_ERROR) { return false; }
          break;
        case 8:
          ontology = *it;
          if (ontology.size() == 0) { return false; }
          break;
        default:
          break;
      }
      ++i;
    }

    ++count;
  }
  in.close();


  if (count < 5) {

    return false;

  } else {

    return true;

  }

}

