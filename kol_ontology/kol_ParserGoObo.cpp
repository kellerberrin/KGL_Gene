//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_ParserGoObo.h"
#include "contrib/kol_GoGraphImpl.h"

#include "kol_GoEnums.h"

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <fstream>

#include <boost/algorithm/string.hpp>


namespace kol = kellerberrin::ontology;
namespace kel = kellerberrin;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Namespace object for the parser implementation.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParserGoOboImpl {

public:

  ParserGoOboImpl() = delete;
  ~ParserGoOboImpl() = delete;

  static kol::GoTermMap parseGoTermFile(const std::string &file_name, const kol::PolicyRelationship& relationship_policy);
  static kol::GoTermMap filterTermMap(const kol::GoTermMap& unfiltered_term_map, const kol::PolicyRelationship& relationship_policy);

private:

  const static constexpr char* SEPARATOR_STR = ": ";  // Separates the "info_id: from the info content."
  const static constexpr char SEPARATOR_SPACE = ' ';
  const static constexpr char* TERM_TOKEN = "[Term]";  // Token that begins a go term.
  const static constexpr char* ID_TOKEN = "id";
  const static constexpr char* ALT_ID_TOKEN = "alt_id";
  const static constexpr char* NAME_TOKEN = "name";
  const static constexpr char* DEF_TOKEN = "def";
  const static constexpr char* NAMESPACE_TOKEN = "namespace";
  const static constexpr char* OBSOLETE_TOKEN = "is_obsolete";
  const static constexpr char* CHILD_TOKEN = "is_a";
  const static constexpr char* RELATIONSHIP_TOKEN = "relationship";

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<kol::GoGraph> kol::ParserGoObo::parseGoFile(const std::string &filename) const {


  auto term_map = ParserGoOboImpl::parseGoTermFile(filename, relationship_policy_);
  auto filtered_term_map = ParserGoOboImpl::filterTermMap(term_map, relationship_policy_);
  std::unique_ptr<GoGraphImpl> graph_impl_ptr(std::make_unique<GoGraphImpl>(filtered_term_map));

  // Return the encapsulated graph.
  return std::make_shared<GoGraph>(std::move(graph_impl_ptr));

}



kol::GoTermMap ParserGoOboImpl::parseGoTermFile(const std::string &file_name, const kol::PolicyRelationship& relationship_policy) {


  std::shared_ptr<kol::GoTermRecord> term_record_ptr(std::make_shared<kol::GoTermRecord>());
  kol::GoTermMap go_term_map;
  enum class GoParserState { FIND_TERM, PROCESS_TERM};

  // Check the relationship policy
  if (not relationship_policy.validPolicy()) {

    kel::ExecEnv::log().error("ParserGoObo::parseGoFile; invalid relationship policy, go file: {}", file_name);
    return go_term_map;

  }

  std::ifstream term_file(file_name);
  if (not term_file.good()) {

    kel::ExecEnv::log().error("ParserGoObo::parseGoFile; problem opening go file: {}", file_name);
    return go_term_map;

  }

  std::string line;
  std::string trimmed_line;
  size_t line_count{0};

  GoParserState parser_state = GoParserState::FIND_TERM;
//  GoTermRecord term_record;
  while (not std::getline(term_file, line).eof()) {

    ++line_count;
    trimmed_line = kel::Utility::trimEndWhiteSpace(line);

    if (parser_state == GoParserState::FIND_TERM) {

      if (trimmed_line == TERM_TOKEN) {

        parser_state = GoParserState::PROCESS_TERM;
        term_record_ptr->clearRecord();

      }

      continue;

    }

    if (trimmed_line.empty())
    {

      if (not term_record_ptr->validRecord()) {

        kel::ExecEnv::log().error("ParserGoObo::parseGoFile; invalid term record at term: {}, line:{}, ",  term_record_ptr->termId(), line_count);
        kel::ExecEnv::log().error("ParserGoObo::parseGoFile; id: {}, name: {}, def: {}, ontology: {}",
                             term_record_ptr->termId(), term_record_ptr->name(), term_record_ptr->definition(), kol::GO::ontologyToString(term_record_ptr->ontology()));

      } else {

        auto [insert_iter, result] = go_term_map.try_emplace(term_record_ptr->termId(), term_record_ptr);
        if (not result) {

          kel::ExecEnv::log().error("ParserGoObo::parseGoFile; cannot insert term, duplicate term id: {}, line:{}",  term_record_ptr->termId(), line_count);

        }

      }
      term_record_ptr = std::make_shared<kol::GoTermRecord>();
      parser_state = GoParserState::FIND_TERM;
      continue;

    }

    auto pos = trimmed_line.find_first_of(SEPARATOR_STR);
    if (pos == std::string::npos) {

      kel::ExecEnv::log().error("ParserGoObo::parseGoFile; expected 'key: value' format, line text : {}, line number: {}",  trimmed_line, line_count);
      parser_state = GoParserState::FIND_TERM;
      continue;

    }

    std::string key = trimmed_line.substr(0, pos);
    std::string value = trimmed_line.substr((pos + 2));

    if (key == OBSOLETE_TOKEN) {

      parser_state = GoParserState::FIND_TERM;
      continue;

    } else if (key == NAMESPACE_TOKEN) {

      if (value == kol::GO::ONTOLOGY_BIOLOGICAL_PROCESS_TEXT) {

        term_record_ptr->ontology(kol::GO::Ontology::BIOLOGICAL_PROCESS);

      } else if (value == kol::GO::ONTOLOGY_MOLECULAR_FUNCTION_TEXT) {

        term_record_ptr->ontology(kol::GO::Ontology::MOLECULAR_FUNCTION);

      } else if (value == kol::GO::ONTOLOGY_CELLULAR_COMPONENT_TEXT) {

        term_record_ptr->ontology(kol::GO::Ontology::CELLULAR_COMPONENT);

      } else {

        kel::ExecEnv::log().error("ParserGoObo::parseGoFile; unexpected namespace: {}",  value);
        parser_state = GoParserState::FIND_TERM;
        continue;

      }

    } else if (key == ID_TOKEN) {

      term_record_ptr->termId(kel::Utility::trimEndWhiteSpace(value));

    } else if (key == ALT_ID_TOKEN) {

      term_record_ptr->altId(kel::Utility::trimEndWhiteSpace(value));

    } else if (key == NAME_TOKEN) {

      term_record_ptr->name(value);

    } else if (key == DEF_TOKEN) {

      term_record_ptr->definition(value);

    } else if (key == CHILD_TOKEN) {

      auto view_vector = kel::Utility::view_tokenizer(value, SEPARATOR_SPACE);
      std::pair<std::string, kol::GO::Relationship> relation{ std::string(view_vector[0]),  kol::GO::Relationship::IS_A};
      term_record_ptr->relations(relation);

    } else if (key == RELATIONSHIP_TOKEN) {

      auto view_vector = kel::Utility::view_tokenizer(value, SEPARATOR_SPACE);
      if (view_vector.size() < 2) {

        kel::ExecEnv::log().error("ParserGoObo::parseGoFile; unable to parse relationship: {}, line: {}", value, line_count);
        parser_state = GoParserState::FIND_TERM;
        continue;

      }

      kol::GO::Relationship relationship = kol::GO::relationshipStringToCode(std::string(view_vector[0]));
      if (relationship == kol::GO::Relationship::REL_ERROR) {

        kel::ExecEnv::log().error("ParserGoObo::parseGoFile; unable to parse relationship: {}, line: {}", value, line_count);
        parser_state = GoParserState::FIND_TERM;
        continue;

      }

      std::pair<std::string, kol::GO::Relationship> relation{ std::string(view_vector[1]),  relationship};
      term_record_ptr->relations(relation);

    } else {

      std::pair<std::string, std::string> attribute{key, value};
      term_record_ptr->attributes(attribute);
    // Additional term informatiomn.

    }

  } // Term record loop.


  return go_term_map;

}


kol::GoTermMap ParserGoOboImpl::filterTermMap(const kol::GoTermMap& unfiltered_term_map, const kol::PolicyRelationship& relationship_policy) {

  kol::GoTermMap filtered_term_map;

  for (auto const& [term_id, term_record_ptr] : unfiltered_term_map) {

    std::vector<std::pair<std::string, kol::GO::Relationship>> filtered_relationships;
    for (auto const& [rel_term, rel_type] : term_record_ptr->relations()) {

      if (relationship_policy.isAllowed(rel_type)) {

        filtered_relationships.emplace_back(rel_term, rel_type);

      }

    }

    if (not filtered_relationships.empty() or kol::GO::isRootTerm(term_id)) {

      std::shared_ptr<kol::GoTermRecord> filtered_record_ptr(std::make_shared<kol::GoTermRecord>(*term_record_ptr));
      filtered_record_ptr->relations(std::move(filtered_relationships));
      auto [insert_iter, result] = filtered_term_map.try_emplace(term_id, filtered_record_ptr);
      if (not result) {

        kel::ExecEnv::log().error("ParserGoObo::filterTermMap; (Duplicate) Unable to insert filtered term: {}", term_id);

      }

    }

  }

  return filtered_term_map;

}


//! Method to parse the go file, should be an OBO file
/*!
  This method will read a Gene Ontology OBO file and add only those relationship
   which are specified to the graph.

*/
std::shared_ptr<kol::GoGraph> kol::ParserGoObo::parseGoFile2(const std::string &filename) const {

  //graph_impl_ptr object to be returned
  std::unique_ptr<GoGraphImpl> graph_impl_ptr(std::make_unique<GoGraphImpl>());

  // Check the relationship policy
  if (not relationship_policy_.validPolicy()) {

    // Return the encapsulated graph_impl_ptr
    return std::make_shared<GoGraph>(std::move(graph_impl_ptr));

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
      //add the term to the graph_impl_ptr term
      graph_impl_ptr->insertTerm(term, name, description, ontology);

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
        graph_impl_ptr->insertTerm(relatedTerm, "name", "description", "ontology");

        //insert edge
        graph_impl_ptr->insertRelationship(term, relatedTerm, relationship);

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

  //call to initialize the graph_impl_ptr's vertex to index maps
  graph_impl_ptr->initMaps();

  // Return the encapsulated graph_impl_ptr;
  return std::make_shared<GoGraph>(std::move(graph_impl_ptr));

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
