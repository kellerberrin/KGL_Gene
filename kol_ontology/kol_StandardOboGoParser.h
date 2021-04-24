/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_STANDARD_OBO_GO_PARSER
#define KGL_STANDARD_OBO_GO_PARSER

#include "kol_GoParserInterface.h"
#include "kol_GoEnums.h"
#include "kol_RelationshipPolicyInterface.h"
#include "kol_StandardRelationshipPolicy.h"
#include "kol_AllowedRelationshipOboGoParser.h"

#include <vector>
#include <string>

#include <boost/tokenizer.hpp>

namespace kellerberrin::ontology {


/*! \class StandardOboGoParser
	\brief A class to parse only is_a or part_of relationships

	 Implements GoParserInterface

*/
class StandardOboGoParser : public AllowedRelationshipOboGoParser {

public:

  StandardOboGoParser() : AllowedRelationshipOboGoParser(StandardRelationshipPolicy()) {}

  ~StandardOboGoParser() override = default;


  //! a method to create a new instance of this class for use in a factory
  /*!
    Creates a new pointer to the parser, used by the factory for go parsers.
  */
  [[nodiscard]] std::unique_ptr<GoParserInterface> clone() const override { return std::make_unique<StandardOboGoParser>(); }


};

} // namespace

#endif
