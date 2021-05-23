/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_STANDARD_OBO_GO_PARSER
#define KGL_STANDARD_OBO_GO_PARSER

#include "kol_ParserGoInterface.h"
#include "kol_GoEnums.h"
#include "kol_PolicyAllowedRelationship.h"
#include "kol_ParserGoAllowedRelationshipObo.h"

#include <vector>
#include <string>

#include <boost/tokenizer.hpp>

namespace kellerberrin::ontology {


/*! \class ParserGoStandardObo
	\brief A class to parse only is_a or part_of relationships

	 Implements ParserGoInterface

*/
class ParserGoStandardObo : public ParserGoAllowedRelationshipObo {

public:

  ParserGoStandardObo() : ParserGoAllowedRelationshipObo(PolicyAllowedRelationship()) {}

  ~ParserGoStandardObo() override = default;


};

} // namespace

#endif
