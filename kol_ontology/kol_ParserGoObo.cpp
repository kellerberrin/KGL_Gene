//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_ParserGoObo.h"

#include "kol_GoEnums.h"

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <iostream>
#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


namespace kol = kellerberrin::ontology;

struct TermRecord {

  std::string term_id;
  std::string name;
  std::string definition;
  kol::GO::Ontology ontology;
  std::vector<std::pair<std::string, kol::GO::Relationship>> relations;
  std::vector<std::string> alt_id;

  void clearRecord();
  bool validRecord() const;

};

void TermRecord::clearRecord() {

  term_id.clear();
  name.clear();
  definition.clear();
  ontology = kol::GO::Ontology::ONTO_ERROR;
  relations.clear();
  alt_id.clear();
}

bool TermRecord::validRecord() const {

  bool valid;
  valid = not term_id.empty();
  valid = valid and not name.empty();
  valid = valid and not definition.empty();
  valid = valid and ontology != kol::GO::Ontology::ONTO_ERROR;
  return valid;

}


std::unique_ptr<kol::GoGraph> kol::ParserGoObo::parseGoFile1(const std::string &file_name) const {

  enum class GoParserState { FIND_TERM, PROCESS_TERM};

  //graph object to be returned
  std::unique_ptr<GoGraph> graph_ptr(std::make_unique<GoGraph>());

  // Check the relationship policy
  if (not relationship_policy_.validPolicy()) {

    ExecEnv::log().error("ParserGoObo::parseGoFile; invalid relationship policy, go file: {}", file_name);
    return graph_ptr;

  }

  std::ifstream term_file(file_name);
  if (not term_file.good()) {

    ExecEnv::log().error("ParserGoObo::parseGoFile; problem opening go file: {}", file_name);
    return graph_ptr;

  }

  std::string line;
  std::string trimmed_line;
  size_t line_count{0};
  size_t term_count{0};
  size_t obsolete_terms{0};
  size_t term_bp{0};
  size_t term_cc{0};
  size_t term_mf{0};

  GoParserState parser_state = GoParserState::FIND_TERM;
  TermRecord term_record;
  while (not std::getline(term_file, line).eof()) {

    ++line_count;
    trimmed_line = Utility::trimEndWhiteSpace(line);

    if (parser_state == GoParserState::FIND_TERM) {

      if (trimmed_line == TERM_TOKEN) {

        parser_state = GoParserState::PROCESS_TERM;
        term_record.clearRecord();
        ++term_count;

      }

      continue;

    }

    if (trimmed_line.empty())
    {

      if (not term_record.validRecord()) {

        ExecEnv::log().error("ParserGoObo::parseGoFile; invalid term record at term: {}, line:{}, ",  term_record.term_id, line_count);
        ExecEnv::log().error("ParserGoObo::parseGoFile; id: {}, name: {}, def: {}, ontology: {}",
                             term_record.term_id, term_record.name, term_record.definition, GO::ontologyToString(term_record.ontology));

      }
      graph_ptr->insertTerm(term_record.term_id, term_record.name, term_record.definition, GO::ontologyToString(term_record.ontology));
      for (auto const& [related_term, relation] : term_record.relations) {
        //insert related terms, there are just stubs to be overwritten later on
        graph_ptr->insertTerm(related_term, "name", "description", "ontology");

        //insert edge
        graph_ptr->insertRelationship(term_record.term_id, related_term, GO::relationshipToString(relation));

      }
      parser_state = GoParserState::FIND_TERM;
      continue;

    }

    auto pos = trimmed_line.find_first_of(SEPARATOR_STR);
    if (pos == std::string::npos) {

      ExecEnv::log().error("ParserGoObo::parseGoFile; expected 'key: value' format, line text : {}, line number: {}",  trimmed_line, line_count);
      parser_state = GoParserState::FIND_TERM;
      continue;

    }

    std::string key = trimmed_line.substr(0, pos);
    std::string value = trimmed_line.substr((pos + 2));

    if (key == OBSOLETE_TOKEN) {

      ++obsolete_terms;
      parser_state = GoParserState::FIND_TERM;
      continue;

    } else if (key == NAMESPACE_TOKEN) {

      if (value == NAMESPACE_BP) {

        term_record.ontology = GO::Ontology::BIOLOGICAL_PROCESS;
        ++term_bp;

      } else if (value == NAMESPACE_MF) {

        term_record.ontology = GO::Ontology::MOLECULAR_FUNCTION;
        ++term_mf;

      } else if (value == NAMESPACE_CC) {

        term_record.ontology = GO::Ontology::CELLULAR_COMPONENT;
        ++term_cc;

      } else {

        ExecEnv::log().error("ParserGoObo::parseGoFile; unexpected namespace: {}",  value);
        parser_state = GoParserState::FIND_TERM;
        continue;

      }

    } else if (key == ID_TOKEN) {

      term_record.term_id = Utility::trimEndWhiteSpace(value);

    } else if (key == ALT_ID_TOKEN) {

      term_record.alt_id.push_back(Utility::trimEndWhiteSpace(value));

    } else if (key == NAME_TOKEN) {

      term_record.name = value;

    } else if (key == DEF_TOKEN) {

      term_record.definition = value;

    } else if (key == CHILD_TOKEN) {

      auto view_vector = Utility::view_tokenizer(value, SEPARATOR_SPACE);
      term_record.relations.emplace_back(std::string(view_vector[0]), GO::Relationship::IS_A);

    } else if (key == RELATIONSHIP_TOKEN) {

      auto view_vector = Utility::view_tokenizer(value, SEPARATOR_SPACE);
      if (view_vector.size() < 2) {

        ExecEnv::log().error("ParserGoObo::parseGoFile; unable to parse relationship: {}, line: {}", value, line_count);
        parser_state = GoParserState::FIND_TERM;
        continue;

      }

      GO::Relationship relation = GO::relationshipStringToCode(std::string(view_vector[0]));
      if (relation == GO::Relationship::REL_ERROR) {

        ExecEnv::log().error("ParserGoObo::parseGoFile; unable to parse relationship: {}, line: {}", value, line_count);
        parser_state = GoParserState::FIND_TERM;
        continue;

      }

      if (relationship_policy_.isAllowed(relation)) {

        term_record.relations.emplace_back(std::string(view_vector[1]), relation);

      }

    } else {

    // Additional term informatiomn.

    }

  }

  ExecEnv::log().info("Parsed GO file: {}, terms: {}, BP: {}, MF: {}, CC: {}, obsolete: {}, lines: {}",
                      file_name, term_count, term_bp, term_mf, term_cc, obsolete_terms, line_count);

  //call to initialize the graph's vertex to index maps
  graph_ptr->initMaps();

  return graph_ptr;

}


//! Method to parse the go file, should be an OBO file
/*!
  This method will read a Gene Ontology OBO file and add only those relationship
   which are specified to the graph.

*/
std::unique_ptr<kol::GoGraph> kol::ParserGoObo::parseGoFile2(const std::string &filename) const {

  //graph object to be returned
  std::unique_ptr<GoGraph> graph(std::make_unique<GoGraph>());

  // Check the relationship policy
  if (not relationship_policy_.validPolicy()) {

    return graph;

  }

  std::ifstream term_file(filename);
  std::string line;

  //at first line
  std::getline(term_file, line);

  //advance to first term
  while (term_file.good() && line != "[Term]") {

    std::getline(term_file, line);
    boost::trim(line);

  }

  while (term_file.good()) {
    //parse a term

    //temp data
    std::string term, name, description, ontology;
    bool isObsolete{false};

    //vector to hold go accession ids of related terms.
    std::vector<std::string> relatedTerms;
    std::vector<std::string> relatedTermRelationship;

    //std::cout << "parsing attributes" << std::endl;
    while (term_file.good() && line != "") {

      std::getline(term_file, line);

      if (line == "") {

        continue;

      }


      std::string attr, value;
      splitWith(line, ":", attr, value);

      //std::cout << "attr: " << attr << " value: " << value << std::endl;
      if (attr == "id") {
        //go accession id
        term = value;

      } else if (attr == "name") {
        //go natural language name/title
        name = value;

      } else if (attr == "namespace") {
        //go ontology, BIOLOGICAL_PROCESS,MOLECULAR_FUNCTION,CELLULAR_COMPONENT
        ontology = value;

      } else if (attr == "def") {
        //go term definition, long form description
        description = value;

      } else if (attr == "is_a") {
        //is_a relationship
        std::string relatedTerm, desc;

        splitWith(value, "!", relatedTerm, desc);

        //add related term to vector for adding later
        relatedTerms.push_back(relatedTerm);

        //add attr as relationship, attr must = "is_a" here
        relatedTermRelationship.push_back(attr);
      } else if (attr == "is_obsolete") {
        //if obsolete set flag.
        isObsolete = true;

      } else if (attr == "relationship") {
        std::string type, toTerm, desc, temp;

        //std::cout << "relationship" << std::endl;
        //std::cout << value << std::endl;
        //extract relationship
        splitWith(value, " ", type, desc);
        //extract term to which related
        splitWith(desc, "!", toTerm, temp);

        //add all types of relationships
        relatedTerms.push_back(toTerm);

        relatedTermRelationship.push_back(type);
      }
    }

    //post process
    //std::cout << term << std::endl;
    //std::cout << name << std::endl;
    //std::cout << description << std::endl;
    //std::cout << ontology << std::endl;


    if (!isObsolete) {
      //add the term to the graph term
      graph->insertTerm(term, name, description, ontology);

      //loop over all related terms
      for (std::size_t i = 0; i < relatedTerms.size(); ++i) {

        //create a temp variable for readability
        std::string relatedTerm = relatedTerms.at(i);
        //std::cout << relatedTerm << "\t";

        //create a temp variable for readability
        std::string relationship = relatedTermRelationship.at(i);
        //std::cout << relationship << std::endl;

        GO::Relationship relationshipType = GO::relationshipStringToCode(relationship);
        if (not relationship_policy_.isAllowed(relationshipType)) {

          continue;

        }

        //insert related terms, there are just stubs to be overwritten later on
        graph->insertTerm(relatedTerm, "name", "description", "ontology");

        //insert edge
        graph->insertRelationship(term, relatedTerm, relationship);

      }//end for, each related term
    }

    //std::cin.get();
    //advance to next term
    while (term_file.good() && line != "[Term]") {
      std::getline(term_file, line);
      boost::trim(line);
    }
    //std::cin.get();
  }

  //call to initialize the graph's vertex to index maps
  graph->initMaps();

  //return the graph pointer
  return graph;
}

//! A method to test if a file fits the accepted format
/*!
  Returns true if the file matches accepted format, false otherwise
*/
bool kol::ParserGoObo::isFileGood(const std::string &filename) const {

  std::ifstream in(filename.c_str());
  if (!in.good()) {
    return false;
  }

  std::string line;
  std::getline(in, line);

  while (in.good() && line != "[Term]") {
    std::getline(in, line);
    boost::trim(line);
  }

  if (!in.good()) {
    return false;
  }

  //parse 3 terms
  std::size_t count = 0;
  while (in.good() && count < 3) {

    while (in.good() && line != "") {
      std::getline(in, line);
      if (line == "") { continue; } // skip any empty lines at the end of input

      std::string attr, value;
      splitWith(line, ":", attr, value);
      if (value == "" || value == "") {
        return false; //if any attributes are empyt return false.
      }
    }

    //Next term
    while (in.good() && line != "[Term]") {
      std::getline(in, line);
      boost::trim(line);
    }
    count += 1;
  }

  if (count < 3 || !in.good()) {
    return false;
  } else {
    return true;
  }

}


//! a helper method
/*!
  splits strings on the given string pattern, splitStr.
*/
void kol::ParserGoObo::splitWith(const std::string &instr,
                                 const std::string &splitStr,
                                 std::string &attr,
                                 std::string &value) const {

  size_t div = instr.find(splitStr);

  if (div != std::string::npos) {

    attr = instr.substr(0, div);
    value = instr.substr(div + 1);
    boost::trim(attr);
    boost::trim(value);

  }

}
