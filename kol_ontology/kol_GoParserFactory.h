/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_GO_PARSER_FACTORY
#define KGL_GO_PARSER_FACTORY

#include "kol_StandardOboGoParser.h"
#include "kol_StandardXmlGoParser.h"
#include "kol_RapidXmlGoParser.h"


#include <vector>
#include <string>

namespace kellerberrin::ontology {

/*! \class GoParserFactory
	\brief A class to return an instance of GoParserInterface at runtime based on an argument.

    STANDARD parsers only parse GO::IS_A and GO::PART_OF go relationships and cannot be modified.
    ALLOWED parsers are initialised to parse GO::IS_A and GO::PART_OF go relationships but this can be modified at runtime.
    XML_RAPID_PARSER parses all relationships.
*/
enum class GoParserType {
  OBO_GO_STANDARD, OBO_GO_ALLOWED, XML_GO_STANDARD, XML_GO_ALLOWED, XML_RAPID_PARSER
};

class GoParserFactory {


public:

  /*!
    This object cannot be created.
  */
  GoParserFactory() = delete;

  ~GoParserFactory() = delete;


  [[nodiscard]] static std::unique_ptr<GoParserInterface> createGoParser(GoParserType parser_type,
                                                                         const RelationshipPolicyInterface &policy = StandardRelationshipPolicy()) {

    switch (parser_type) {

      case GoParserType::OBO_GO_ALLOWED:
        return std::make_unique<AllowedRelationshipOboGoParser>(policy);

      case GoParserType::XML_GO_STANDARD:
        return std::make_unique<StandardXmlGoParser>();

      case GoParserType::XML_GO_ALLOWED:
        return std::make_unique<AllowedRelationshipXmlGoParser>(policy);

      case GoParserType::XML_RAPID_PARSER:
        return std::make_unique<RapidXmlGoParser>();

      default:
      case GoParserType::OBO_GO_STANDARD:
        return std::make_unique<StandardOboGoParser>();

    }

  }


};


} // namespace


#endif