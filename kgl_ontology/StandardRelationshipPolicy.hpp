/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef STANDARD_RELATIONSHIP_POLICY
#define STANDARD_RELATIONSHIP_POLICY

#include <RelationshipPolicyInterface.hpp>
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

		_relationshipMap[GO::IS_A] = true;
		_relationshipMap[GO::PART_OF] = true;

	}
  ~StandardRelationshipPolicy() override = default;

  std::unique_ptr<RelationshipPolicyInterface> clone() const override { return std::make_unique<StandardRelationshipPolicy>(); }

	//! a method to test if a relatinoship is allowed or not
	/*!
		tests if the relationship is allowed. Overridden to fulfill the RelationshipPolicyInterface
	*/
	bool isAllowed(GO::Relationship relationship) const override {return _relationshipMap.find(relationship) != _relationshipMap.end(); }

private:
	//! a map of relationships to bool
    /*! maps a relationship to a bool. Boost unordered map give constant time find. */
	OntologyMapType<GO::Relationship,bool> _relationshipMap;

};
#endif
