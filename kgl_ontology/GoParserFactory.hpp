/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef GO_PARSER_FACTORY
#define GO_PARSER_FACTORY

#include <StandardOboGoParser.hpp>
#include <StandardXmlGoParser.hpp>
#include <RapidXmlGoParser.hpp>


#include <vector>
#include <string>


/*! \class GoParserFactory
	\brief A class to return an instance of GoParserInterface at runtime based on an argument.

    STANDARD parsers only parse GO::IS_A and GO::PART_OF go relationships and cannot be modified.
    ALLOWED parsers are initialised to parse GO::IS_A and GO::PART_OF go relationships but this can be modified at runtime.
    RAPID_XML_PARSER parses all relationships.
*/
enum class GoParserType { OBO_GO_STANDARD, OBO_GO_ALLOWED, XML_GO_STANDARD, XML_GO_ALLOWED, RAPID_XML_PARSER };
class GoParserFactory {


public:

	/*!
		This object cannot be created.
	*/
	GoParserFactory() = delete;
  ~GoParserFactory() = delete;

	//! A Method to add a parser to the factory.
	/*!
		This method adds a pointer to a parser and a string to the factory.
		  This string is used to query the appropriate parser.
	*/

  static std::unique_ptr<GoParserInterface> createGoParser(GoParserType parser_type) {

    switch(parser_type) {

      case GoParserType::OBO_GO_STANDARD:
        return std::make_unique<StandardOboGoParser>();

      case GoParserType::OBO_GO_ALLOWED:
        return std::make_unique<AllowedRelationshipOboGoParser>(StandardRelationshipPolicy());

      case GoParserType::XML_GO_STANDARD:
        return std::make_unique<StandardXmlGoParser>();

      case GoParserType::XML_GO_ALLOWED:
        return std::make_unique<AllowedRelationshipXmlGoParser>(StandardRelationshipPolicy());

      case GoParserType::RAPID_XML_PARSER:
        return std::make_unique<RapidXmlGoParser>();

    }

  }


};
#endif