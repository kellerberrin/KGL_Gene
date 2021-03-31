/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef ALLOWED_SET_RELATIONSHIP_POLICY
#define ALLOWED_SET_RELATIONSHIP_POLICY

#include <vector>
#include <GoEnums.hpp>
#include <RelationshipPolicyInterface.hpp>
#include <boost/unordered_map.hpp>

/*! \class AllowedSetRelationshipPolicy
	\brief A class to allow only a set of relationships

	A class to allow only certain relationships in the go graph. It uses a set of 
	enums to restric the types of relationships considered in a graph.

*/
class AllowedSetRelationshipPolicy: public RelationshipPolicyInterface{

public:
	
	//! A constructor
	/*!
		Creates the default(empty) AllowedSetRelationshipPolicy
	*/
	AllowedSetRelationshipPolicy() = default;
  AllowedSetRelationshipPolicy(const AllowedSetRelationshipPolicy& copy) = default;
  ~AllowedSetRelationshipPolicy() override = default;

	//! A parameterized constructor
	/*!
		Creats the AllowedSetRelationshipPolicy using a list(vector) of relationships to allow
	*/
	explicit AllowedSetRelationshipPolicy(const std::vector<GO::Relationship>& relationships) {

		for(auto const& relation : relationships) {

			_relationshipMap[relation] = true;

		}

	}

	std::unique_ptr<RelationshipPolicyInterface> clone() const override { return std::make_unique<AllowedSetRelationshipPolicy>(*this); }

	//! a method to test if a relatinoship is allowed or not
	/*!
		tests if the relationship is allowed. Overridden to fulfill the RelationshipPolicyInterface
	*/
	[[nodiscard]] bool isAllowed(GO::Relationship relationship) const override {

	  return _relationshipMap.find(relationship) != _relationshipMap.end();

	}


	//! a method to add a relationship to the set of relationships allowed
	/*!
		adds a relationship to the set of relationships allowed by setting its mapped value to true
	*/
	void addRelationship(GO::Relationship relationship) { _relationshipMap[relationship] = true; }

	//! a method to add a relationship to the set of relationships allowed
	/*!
	adds a relationship to the set of relationships allowed by setting its mapped value to true
	*/
	void addRelationship(const std::string &relString) {

		GO::Relationship relationship = GO::relationshipStringToCode(relString);

		if (relationship != GO::Relationship::REL_ERROR){

			_relationshipMap[relationship] = true;

		}

	}


	//! a method to determine if the Policy is empty
	/*!
		Determines if the Policy is empty
	*/
	inline bool isEmpty() { return _relationshipMap.empty(); }

private:
	//! a map of relationships to bool
    /*! maps a relationship to a bool. Boost unordered map give constant time find. */
	OntologyMapType<GO::Relationship,bool> _relationshipMap;

};
#endif
