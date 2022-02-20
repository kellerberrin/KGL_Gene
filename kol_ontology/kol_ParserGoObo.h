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
  [[nodiscard]] std::shared_ptr<GoGraph> parseGoFile(const std::string &filename) const override;

  //! A method to test if a file fits the accepted format
  /*!
    Returns true if the file matches accepted format, false otherwise
  */

  [[nodiscard]] bool isFileGood(const std::string &filename) const override;

private:


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

  std::shared_ptr<GoGraph>  parseGoFile2(const std::string &filename) const;

};

} // namespace

#endif
