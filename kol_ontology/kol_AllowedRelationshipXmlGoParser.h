/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ALLOWED_RELATIONSHIP_XML_GO_PARSER
#define KGL_ALLOWED_RELATIONSHIP_XML_GO_PARSER

#include "kol_GoParserInterface.h"
#include "kol_GoEnums.h"
#include "kol_RelationshipPolicyInterface.h"

#include <vector>
#include <string>


namespace kellerberrin::ontology {


/*! \class AllowedRelationshipXmlGoParser
	\brief A class to parse only a specified set of relationships

	This class will read a Gene Ontology XML file and add only those relationship
	 which are specified to the graph. The most important method of this class if the
	 parseGoFile which takes the file name as a parameter.

	 Implements GoParserInterface

*/
class AllowedRelationshipXmlGoParser : public GoParserInterface {

public:

  //! A parameterized constructor
  /*!
    constructor that sets the policy
  */
  explicit AllowedRelationshipXmlGoParser(const RelationshipPolicyInterface &policy)
      : relationship_policy_ptr_(policy.clone()) {}

  AllowedRelationshipXmlGoParser() = delete; // Must provide a policy.
  ~AllowedRelationshipXmlGoParser() override = default;

  //! a method to create a new instance of this class for use in a factory
  /*!
    creats a new pointer to the parser, used by the factory for go parsers.
  */
  [[nodiscard]] std::unique_ptr<GoParserInterface> clone() const override {

    return std::make_unique<AllowedRelationshipXmlGoParser>(*relationship_policy_ptr_);

  }//end method clone


  //! Method to parse the go file, should be an XML file
  /*!
    This method will read a Gene Ontology XML file and add only those relationship
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

};

} // namespace

#endif
