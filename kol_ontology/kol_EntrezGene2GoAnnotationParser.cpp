//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_EntrezGene2GoAnnotationParser.h"

#include <iostream>
#include <boost/tokenizer.hpp>

namespace kol = kellerberrin::ontology;



//! An interface method for parsing an annotation file.
/*!
  This method takes a filename as in put and returns a pointer to an
    AnnotationData object. This method fulfills part of the interface contract.
*/
std::unique_ptr<kol::AnnotationData> kol::EntrezGene2GoAnnotationParser::parseAnnotationFile(const std::string &filename) const  {

  std::unique_ptr<AnnotationData> annoData(std::make_unique<AnnotationData>());

  // Check that the supplied policy is valid.
  if (not policy_ptr_->isValid()) {

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
    getline(in, line);
    if (line[0] == '#') { continue; }

    //split line
    tokenizer tokens(line, tab_sep);
    //set iterator to first token
    it = tokens.begin();

    //always check if empty (must have in linux)
    if (it == tokens.end()) { continue; }
    //iterator at 0

    //temp storage variables
    std::string taxon, geneName, goTerm, evidenceCode;

    //iterator at 0
    //set taxon
    taxon = *it;

    //only use annotations for the specified taxonomy
    //if(taxon.compare(taxonomy) != 0){continue;}

    std::advance(it, 1);
    //iterator at 1
    //set gene id
    geneName = *it;

    std::advance(it, 1);
    //iterator at 2
    //set go term
    goTerm = *it;

    std::advance(it, 1);
    //iterator at 3
    //set evidence code
    evidenceCode = *it;

    //add gene to go association to the database, if the evidence is allowed
    if (policy_ptr_->isAllowed(GO::evidenceStringToCode(evidenceCode))) {
      //add gene to go association to the database
      annoData->addAssociation(geneName, goTerm, evidenceCode);
    }
  }
  in.close();

  return annoData;
}


//! A method for checking if a file exists and is formatted correctly.
/*!
  This function checks that the file exists and its format can be recognized.
*/
bool kol::EntrezGene2GoAnnotationParser::isFileGood(const std::string &fileName) const {

  std::ifstream in(fileName.c_str());
  if (!in.good()) {
    return false;
  }

  //Tokenizer type
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> tab_sep("\t", "", boost::keep_empty_tokens);
  tokenizer::iterator it;
  std::string line;

  std::size_t count = 0;

  //main loop, each line
  while (in.good() && count < 5) {
    //get next line in 'in' file stream
    getline(in, line);
    if (line[0] == '#' or line[0] == '!') { continue; }

    //split line
    tokenizer tokens(line, tab_sep);
    //set iterator to first token
    it = tokens.begin();

    //always check if empty (must have in linux)
    if (it == tokens.end()) { continue; }
    //iterator at 0

    //temp storage variables
    std::string taxon, geneName, goTerm, evidenceCode;
    taxon = *it;
    if (taxon.size() == 0) { return false; }

    std::advance(it, 1);
    geneName = *it;
    if (geneName.size() == 0) { return false; }

    std::advance(it, 1);
    goTerm = *it;
    if (goTerm.size() == 0) { return false; }           // disallow empty go
    if (goTerm.substr(0, 3) != "GO:") { return false; } // disallow bad go term

    std::advance(it, 1);
    evidenceCode = *it;
    if (evidenceCode.size() == 0) { return false; }
    if (GO::evidenceStringToCode(evidenceCode) == GO::EvidenceCode::ECODE_ERROR) { return false; }

    ++count;
  }
  in.close();


  if (count < 5) {
    return false;
  } else {
    return true;
  }
}

