/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_ALLOWED_RELATIONSHIP_OBO_GO_PARSER
#define KOL_ALLOWED_RELATIONSHIP_OBO_GO_PARSER

#include "kol_ParserGoInterface.h"
#include "kol_GoEnums.h"
#include "kol_PolicyRelationship.h"
#include <vector>
#include <string>


namespace kellerberrin::ontology {

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
  [[nodiscard]] std::shared_ptr<GoGraphImpl> parseGoFile(const std::string &filename) const override;

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

  GoTermMap parseGoTermFile(const std::string &file_name) const;
  std::shared_ptr<GoGraphImpl>  parseGoFile2(const std::string &filename) const;
  GoTermMap filterTermMap(const GoTermMap& unfiltered_term_map, const PolicyRelationship& relationship_policy) const;

};

} // namespace

#endif
