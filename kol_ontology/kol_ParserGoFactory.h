/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_GO_PARSER_FACTORY
#define KOL_GO_PARSER_FACTORY

#include "kol_ParserGoObo.h"
#include "kol_ParserGoXml.h"
#include "kol_ParserGoRapidXml.h"


#include <vector>
#include <string>

namespace kellerberrin::ontology {

/*! \class ParserGoFactory
	\brief A class to return an instance of ParserGoInterface at runtime based on an argument.

    STANDARD parsers only parse GO::IS_A and GO::PART_OF go relationships and cannot be modified.
    ALLOWED parsers are initialised to parse GO::IS_A and GO::PART_OF go relationships but this can be modified at runtime.
    PARSER_GO_XML_RAPID parses all relationships.
*/
enum class ParserGoType { PARSER_GO_OBO, PARSER_GO_XML, PARSER_GO_XML_RAPID };

class ParserGoFactory {


public:

  /*!
    This object cannot be created.
  */
  ParserGoFactory() = delete;
  ~ParserGoFactory() = delete;


  [[nodiscard]] static std::unique_ptr<ParserGoInterface> createGoParser(ParserGoType parser_type,
                                                                         const PolicyRelationship &policy = PolicyRelationship()) {

    switch (parser_type) {

      case ParserGoType::PARSER_GO_XML:
        return std::make_unique<ParserGoXml>(policy);

      case ParserGoType::PARSER_GO_XML_RAPID:
        return std::make_unique<GoParserRapidXml>();

      default:
      case ParserGoType::PARSER_GO_OBO:
        return std::make_unique<ParserGoObo>(policy);

    }

  }


};


} // namespace


#endif