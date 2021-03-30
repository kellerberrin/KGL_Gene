/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef STANDARD_OBO_GO_PARSER
#define STANDARD_OBO_GO_PARSER

#include <GoParserInterface.hpp>
#include <GoEnums.hpp>
#include <RelationshipPolicyInterface.hpp>
#include <StandardRelationshipPolicy.hpp>
#include <AllowedRelationshipOboGoParser.hpp>

#include <vector>
#include <string>

#include <boost/tokenizer.hpp>


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

	// No function
  void setPolicy(const RelationshipPolicyInterface&) override {}

};
#endif
