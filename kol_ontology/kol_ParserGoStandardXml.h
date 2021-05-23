/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_STANDARD_XML_GO_PARSER
#define KGL_STANDARD_XML_GO_PARSER

#include "kol_ParserGoInterface.h"
#include "kol_GoEnums.h"
#include "kol_PolicyAllowedRelationship.h"
#include "kol_ParserGoAllowedRelationshipXml.h"

#include <vector>
#include <string>

#include <xml/rapidxml_utils.h>
#include <xml/rapidxml.h>

#include <boost/tokenizer.hpp>


namespace kellerberrin::ontology {

/*! \class GoParserStandardXml
	\brief A class to parse only is_a or part_of relationships

	 Implements ParserGoInterface

*/
class GoParserStandardXml : public ParserGoAllowedRelationshipXml {

public:

  GoParserStandardXml() : ParserGoAllowedRelationshipXml(PolicyAllowedRelationship()) {}

  ~GoParserStandardXml() override = default;


};

} // namespace

#endif
