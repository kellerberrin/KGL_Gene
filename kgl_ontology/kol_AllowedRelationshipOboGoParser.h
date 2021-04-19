/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ALLOWED_RELATIONSHIP_OBO_GO_PARSER
#define KGL_ALLOWED_RELATIONSHIP_OBO_GO_PARSER

#include "kol_GoParserInterface.h"
#include "kol_GoEnums.h"
#include "kol_RelationshipPolicyInterface.h"
#include <vector>
#include <string>


namespace kellerberrin::ontology {


/*! \class AllowedRelationshipOboGoParser
	\brief A class to parse only a specified set of relationships

	This class will read a Gene Ontology OBO file and add only those relationships
	 which are specified to the graph. The most important method of this class if the
	 parseGoFile which takes the file name as a parameter.

	 Implements GoParserInterface

*/
class AllowedRelationshipOboGoParser : public GoParserInterface {

public:

  //! A parameterized constructor
  /*!
    constructor that sets the policy
  */
  explicit AllowedRelationshipOboGoParser(const RelationshipPolicyInterface &policy)
      : relationship_policy_ptr_(policy.clone()) {}

  AllowedRelationshipOboGoParser() = delete; // Must provide a policy.
  ~AllowedRelationshipOboGoParser() override = default;

  //! a method to create a new instance of this class for use in a factory
  /*!
    creats a new pointer to the parser, used by the factory for go parsers.
  */
  [[nodiscard]] std::unique_ptr<GoParserInterface> clone() const override {

    return std::make_unique<AllowedRelationshipOboGoParser>(*relationship_policy_ptr_);

  }


  //! Method to parse the go file, should be an OBO file
  /*!
    This method will read a Gene Ontology OBO file and add only those relationship
     which are specified to the graph.

  */
  [[nodiscard]] std::unique_ptr<GoGraph> parseGoFile(const std::string &filename) const override;

  //! A method to test if a file fits the accepted format
  /*!
    Returns true if the file matches accepted format, false otherwise
  */

  [[nodiscard]] bool isFileGood(const std::string &filename) const override;

private:

  //! A RelationshipPolicyInterface
  /*! This RelationshipPolicyInterface holds the relationships to be allowed during parsing */
  std::unique_ptr<RelationshipPolicyInterface> relationship_policy_ptr_;


  //! a helper method
  /*!
    splits strings on the given string pattern, splitStr.
  */
  void splitWith(const std::string &instr,
                 const std::string &splitStr,
                 std::string &attr,
                 std::string &value) const;

};

} // namespace

#endif
