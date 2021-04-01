/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef STANDARD_XML_GO_PARSER
#define STANDARD_XML_GO_PARSER

#include <GoParserInterface.h>
#include <GoEnums.h>
#include <RelationshipPolicyInterface.h>
#include <StandardRelationshipPolicy.h>
#include <AllowedRelationshipXmlGoParser.h>

#include <vector>
#include <string>

#include <xml/rapidxml_utils.h>
#include <xml/rapidxml.h>

#include <boost/tokenizer.hpp>


/*! \class StandardXmlGoParser
	\brief A class to parse only is_a or part_of relationships

	 Implements GoParserInterface

*/
class StandardXmlGoParser : public AllowedRelationshipXmlGoParser {

public:

  StandardXmlGoParser() : AllowedRelationshipXmlGoParser(StandardRelationshipPolicy()) {}
  ~StandardXmlGoParser() override = default;

  //! a method to create a new instance of this class for use in a factory
  /*!
    Creates a new pointer to the parser, used by the factory for go parsers.
  */

	[[nodiscard]] std::unique_ptr<GoParserInterface> clone() const override { return std::make_unique<StandardXmlGoParser>(); }


};
#endif
