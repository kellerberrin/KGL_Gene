/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ALLOWED_RELATIONSHIP_OBO_GO_PARSER
#define KGL_ALLOWED_RELATIONSHIP_OBO_GO_PARSER

#include "kol_ParserGoInterface.h"
#include "kol_GoEnums.h"
#include "kol_PolicyRelationship.h"
#include <vector>
#include <string>


namespace kellerberrin::ontology {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GoTermRecord {

public:

  GoTermRecord() = default;
  GoTermRecord(const GoTermRecord&) = default;
  ~GoTermRecord() = default;

  // Getters
  [[nodiscard]] const std::string& termId() const { return term_id_; }
  [[nodiscard]] const std::string& name() const { return name_; }
  [[nodiscard]] const std::string& definition() const { return definition_; }
  [[nodiscard]] GO::Ontology ontology() const { return ontology_; }
  [[nodiscard]] const std::vector<std::pair<std::string, GO::Relationship>>& relations() const { return relations_; }
  [[nodiscard]] const std::vector<std::string>& altId() const { return alt_id_; }
  [[nodiscard]] const std::vector<std::pair<std::string, std::string>>& attributes() const { return attributes_; }

  // Setters.
  void termId(const std::string& term_id) { term_id_ = term_id; }
  void name(const std::string& name) { name_ = name; }
  void definition(const std::string& definition) { definition_ = definition; }
  void ontology(GO::Ontology ontology) { ontology_ = ontology; }
  void relations(std::pair<std::string, GO::Relationship> relation) { relations_.push_back(std::move(relation)); }
  void altId(const std::string& alt_id) { alt_id_.push_back(alt_id); }
  void attributes(std::pair<std::string, std::string> attribute) { attributes_.push_back(std::move(attribute)); }

  void clearRecord();
  [[nodiscard]] bool validRecord() const;

private:

  std::string term_id_;
  std::string name_;
  std::string definition_;
  GO::Ontology ontology_{GO::Ontology::ONTO_ERROR};
  std::vector<std::pair<std::string, GO::Relationship>> relations_;
  std::vector<std::string> alt_id_;
  std::vector<std::pair<std::string, std::string>> attributes_;  // key:value pairs is all misc attributes.

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*! \class ParserGoObo
	\brief A class to parse only a specified set of relationships

	This class will read a Gene Ontology OBO file and add only those relationships
	 which are specified to the graph. The most important method of this class if the
	 parseGoFile which takes the file name as a parameter.

	 Implements ParserGoInterface

*/
class ParserGoObo : public ParserGoInterface {

public:

  //! A parameterized constructor
  /*!
    constructor that sets the policy
  */
  explicit ParserGoObo(const PolicyRelationship &policy) : relationship_policy_(policy) {}
  ParserGoObo() = delete; // Must provide a policy.
  ~ParserGoObo() override = default;


  //! Method to parse the go file, should be an OBO file
  /*!
    This method will read a Gene Ontology OBO file and add only those relationship
     which are specified to the graph.

  */
  [[nodiscard]] std::shared_ptr<GoGraph> parseGoFile(const std::string &filename) const override {

    parseGoFile2(filename);
    return parseGoFile1(filename);

  }

  //! A method to test if a file fits the accepted format
  /*!
    Returns true if the file matches accepted format, false otherwise
  */

  [[nodiscard]] bool isFileGood(const std::string &filename) const override;

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
  const static constexpr char* NAMESPACE_BP = "biological_process";
  const static constexpr char* NAMESPACE_MF = "molecular_function";
  const static constexpr char* NAMESPACE_CC = "cellular_component";

  //! A PolicyRelationshipInterface
  /*! This PolicyRelationshipInterface holds the relationships to be allowed during parsing */

  PolicyRelationship relationship_policy_;


  //! a helper method
  /*!
    splits strings on the given string pattern, splitStr.
  */
  void splitWith(const std::string &instr,
                 const std::string &splitStr,
                 std::string &attr,
                 std::string &value) const;

  std::shared_ptr<GoGraph> parseGoFile1(const std::string &file_name) const;
  std::shared_ptr<GoGraph>  parseGoFile2(const std::string &filename) const;

};

} // namespace

#endif
