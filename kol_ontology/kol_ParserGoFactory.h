/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_GO_PARSER_FACTORY
#define KGL_GO_PARSER_FACTORY

#include "kol_ParserGoStandardObo.h"
#include "kol_ParserGoStandardXml.h"
#include "kol_ParserGoRapidXml.h"


#include <vector>
#include <string>

namespace kellerberrin::ontology {

/*! \class ParserGoFactory
	\brief A class to return an instance of ParserGoInterface at runtime based on an argument.

    STANDARD parsers only parse GO::IS_A and GO::PART_OF go relationships and cannot be modified.
    ALLOWED parsers are initialised to parse GO::IS_A and GO::PART_OF go relationships but this can be modified at runtime.
    XML_RAPID_PARSER parses all relationships.
*/
enum class ParserGoType { OBO_GO_STANDARD, OBO_GO_ALLOWED, XML_GO_STANDARD, XML_GO_ALLOWED, XML_RAPID_PARSER };

class ParserGoFactory {


public:

  /*!
    This object cannot be created.
  */
  ParserGoFactory() = delete;
  ~ParserGoFactory() = delete;


  [[nodiscard]] static std::unique_ptr<ParserGoInterface> createGoParser(ParserGoType parser_type,
                                                                         const PolicyAllowedRelationship &policy = PolicyAllowedRelationship()) {

    switch (parser_type) {

      case ParserGoType::OBO_GO_ALLOWED:
        return std::make_unique<ParserGoAllowedRelationshipObo>(policy);

      case ParserGoType::XML_GO_STANDARD:
        return std::make_unique<GoParserStandardXml>();

      case ParserGoType::XML_GO_ALLOWED:
        return std::make_unique<ParserGoAllowedRelationshipXml>(policy);

      case ParserGoType::XML_RAPID_PARSER:
        return std::make_unique<GoParserRapidXml>();

      default:
      case ParserGoType::OBO_GO_STANDARD:
        return std::make_unique<ParserGoStandardObo>();

    }

  }


};


} // namespace


#endif