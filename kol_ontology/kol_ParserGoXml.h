/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ALLOWED_RELATIONSHIP_XML_GO_PARSER
#define KGL_ALLOWED_RELATIONSHIP_XML_GO_PARSER

#include "kol_ParserGoInterface.h"
#include "kol_GoEnums.h"
#include "kol_PolicyRelationship.h"

#include <vector>
#include <string>


namespace kellerberrin::ontology {


/*! \class ParserGoXml
	\brief A class to parse only a specified set of relationships

	This class will read a Gene Ontology XML file and add only those relationship
	 which are specified to the graph. The most important method of this class if the
	 parseGoFile which takes the file name as a parameter.

	 Implements ParserGoInterface

*/
class ParserGoXml : public ParserGoInterface {

public:

  //! A parameterized constructor
  /*!
    constructor that sets the policy
  */
  explicit ParserGoXml(const PolicyRelationship &policy) : relationship_policy_(policy) {}

  ParserGoXml() = delete; // Must provide a policy.
  ~ParserGoXml() override = default;


  //! Method to parse the go file, should be an XML file
  /*!
    This method will read a Gene Ontology XML file and add only those relationship
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

};

} // namespace

#endif
