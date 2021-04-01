/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef STANDARD_RELATIONSHIP_POLICY
#define STANDARD_RELATIONSHIP_POLICY

#include <RelationshipPolicyInterface.h>
#include <boost/unordered_map.hpp>

/*! \class StandardRelationshipPolicy
	\brief A class to allow only a set of relationships

	A class to allow only certain relationships in the go graph. It uses a set of 
	enums to restrict the types of relationships considered in a graph.

*/
class StandardRelationshipPolicy: public RelationshipPolicyInterface{

public:
	
	//! A constructor
	/*!
		Creates the default StandardRelationshipPolicy
	*/
	StandardRelationshipPolicy() {

		_relationshipMap[GO::Relationship::IS_A] = true;
		_relationshipMap[GO::Relationship::PART_OF] = true;

	}
  ~StandardRelationshipPolicy() override = default;

  [[nodiscard]] std::unique_ptr<RelationshipPolicyInterface> clone() const override { return std::make_unique<StandardRelationshipPolicy>(); }

	//! a method to test if a relatinoship is allowed or not
	/*!
		tests if the relationship is allowed. Overridden to fulfill the RelationshipPolicyInterface
	*/
	[[nodiscard]] bool isAllowed(GO::Relationship relationship) const override {return _relationshipMap.find(relationship) != _relationshipMap.end(); }

  //! a method to determine if the Policy is valid
  /*!
    Determines if the Policy is valid.
    Must contain both Relationship::IS_A and Relationship::PART_OF.
    Cannot contain Relationship::REL_ERROR.
    Trivially true in the standard case.
  */

  [[nodiscard]] bool validPolicy() const override { return true; }


private:
	//! a map of relationships to bool
    /*! maps a relationship to a bool. Boost unordered map give constant time find. */
	OntologyMapType<GO::Relationship,bool> _relationshipMap;

};
#endif
