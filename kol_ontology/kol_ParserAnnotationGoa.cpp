//
// Created by kellerberrin on 26/5/21.
//
#include "kol_ParserAnnotationGoa.h"
#include "kol_NewAnnotationData.h"

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <iostream>
#include <boost/tokenizer.hpp>


namespace kol = kellerberrin::ontology;


std::shared_ptr<kol::AnnotationData> kol::ParserAnnotationGoa::parseAnnotationFile(const std::string &file_name) const {

  return parseAnnotationFileNew(file_name);
//  return parseAnnotationFileOld(file_name);

}

std::shared_ptr<kol::AnnotationData> kol::ParserAnnotationGoa::parseAnnotationFileNew(const std::string &file_name) const {

  std::shared_ptr<AnnotationData> annoData(std::make_unique<AnnotationData>());
  std::shared_ptr<AnnotationDataNew> annoData1(std::make_unique<AnnotationDataNew>());


  // Check that the supplied policy is valid.
  if (not _policy.isValid()) {

    return annoData;

  }

  std::ifstream anno_file(file_name);
  if (not anno_file.good()) {

    ExecEnv::log().error("ParserAnnotationGoa::parseAnnotationFile; problem opening annotation file: {}", file_name);
    return annoData;

  }

  std::string line;
  size_t line_count{0};
  size_t record_count{0};
  size_t annotation_count{0};
  size_t count_bf{0};
  size_t count_mf{0};
  size_t count_cc{0};
  size_t count_test{0};
  size_t count_not{0};
  const std::string term("GO:0002269");
  const std::string parent_term("GO:0002035");

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

    ++record_count;
    std::string evidence_code_text(view_vector[FIELD_OFFSET_EVIDENCE_CODE]);
    GO::EvidenceCode evidence_code = GO::evidenceStringToCode(evidence_code_text);
    if (_policy.isAllowed(evidence_code)) {
      //add gene to go association to the database
      if (not view_vector[FIELD_OFFSET_QUALIFIER].empty()) {

        auto qualifier = Utility::toupper(std::string(view_vector[FIELD_OFFSET_QUALIFIER]));
        if (qualifier.find(NOT_QUALIFIER) != std::string::npos) {

          ++count_not;
          continue;

        }

      }
      ++annotation_count;
      std::string uniprot_id(view_vector[FIELD_OFFSET_UNIPROT_ID]);
      std::string term_id(view_vector[FIELD_OFFSET_GO_ID]);
      std::string ontology_code(view_vector[FIELD_OFFSET_ONTOLOGY]);
      GO::Ontology go_ontology;

      if (ontology_code == ONTOLOGY_BP_CODE) {
        ++count_bf;
        go_ontology = GO::Ontology::BIOLOGICAL_PROCESS;
      } else if (ontology_code == ONTOLOGY_MF_CODE) {
        ++count_mf;
        go_ontology = GO::Ontology::MOLECULAR_FUNCTION;
      } else if (ontology_code == ONTOLOGY_CC_CODE) {
        ++count_cc;
        go_ontology = GO::Ontology::CELLULAR_COMPONENT;
      } else {
        ExecEnv::log().error("ParserAnnotationGoa::parseAnnotationFile; line: {}, bad ontology code: {}", line_count, ontology_code);
        continue;
      }
      if (term_id == term or term_id == parent_term) {
        ++count_test;
      }

      annoData->addAssociation(uniprot_id, term_id, evidence_code_text);
      annoData1->addAssociation(uniprot_id, term_id, go_ontology, evidence_code);

    }

  }

  ExecEnv::log().info("AnnotationFile: {}, parsed records: {}, annotations: {}, BP: {}, MF: {}, CC: {}, Test: {}, NOT: {}",
                      file_name, record_count, annotation_count, count_bf, count_mf, count_cc, count_test, count_not);

  return annoData;

}

//! An interface method for parsing an annotation file.
/*!
  This method takes a filename as in put and returns a pointer to an
    AnnotationData object. This method fulfills part of the interface contract.
*/
std::shared_ptr<kol::AnnotationData> kol::ParserAnnotationGoa::parseAnnotationFileOld(const std::string &filename) const {

  std::shared_ptr<AnnotationData> annoData(std::make_unique<AnnotationData>());

  // Check that the supplied policy is valid.
  if (not _policy.isValid()) {

    return annoData;

  }

  //open afile stream
  std::ifstream in(filename.c_str());

  //Tokenizer type
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  //declares a tab separator variable
  boost::char_separator<char> tab_sep("\t", "", boost::keep_empty_tokens);

  //An iterator for the tokens
  tokenizer::iterator it;

  //string variable for each line
  std::string line;

  //main loop, each line
  while (in.good()) {
    //get next line in 'in' file stream
    std::getline(in, line);

    //split line
    tokenizer tokens(line, tab_sep);
    //set iterator to first token
    it = tokens.begin();

    //always check if empty (must have in linux)
    if (it == tokens.end()) { continue; }
    //iterator at 0

    //skip comments if line is not empty and starts with !
    if (line.at(0) == '!') { continue; }


    std::string database, geneName, qualifierStr, goStr, evidenceCode, ontology;

    //database,first field
    database = *it;

    std::advance(it, 1);
    //iterator at 1
    geneName = *it;


    std::advance(it, 2);
    //iterator at 3
    qualifierStr = *it;


    std::advance(it, 1);
    //iterator at 4
    goStr = *it;


    std::advance(it, 2);
    //iterator at 6
    evidenceCode = *it;


    std::advance(it, 2);
    //iterator at 8
    ontology = *it;


    if (_policy.isAllowed(GO::evidenceStringToCode(evidenceCode))) {
      //add gene to go association to the database
      annoData->addAssociation(geneName, goStr, evidenceCode);

    }

  }//end while, each line
  in.close();

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

