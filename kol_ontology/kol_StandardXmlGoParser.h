/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_STANDARD_XML_GO_PARSER
#define KGL_STANDARD_XML_GO_PARSER

#include "kol_GoParserInterface.h"
#include "kol_GoEnums.h"
#include "kol_RelationshipPolicyInterface.h"
#include "kol_StandardRelationshipPolicy.h"
#include "kol_AllowedRelationshipXmlGoParser.h"

#include <vector>
#include <string>

#include <xml/rapidxml_utils.h>
#include <xml/rapidxml.h>

#include <boost/tokenizer.hpp>


namespace kellerberrin::ontology {

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

} // namespace

#endif
