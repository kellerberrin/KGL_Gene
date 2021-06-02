/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_ALLOWED_SET_RELATIONSHIP_POLICY
#define KOL_ALLOWED_SET_RELATIONSHIP_POLICY

#include <vector>
#include "kol_GoEnums.h"


namespace kellerberrin::ontology {

/*! \class PolicyRelationship
	\brief A class to allow only a set of relationships

	A class to allow only certain relationships in the go graph. It uses a set of 
	enums to restric the types of relationships considered in a graph.

*/
class PolicyRelationship  {

public:

  //! A constructor
  /*!
    Creates the default PolicyRelationship
  */
  PolicyRelationship() { addRelationships(GO::standardRelationships()); }
  ~PolicyRelationship()  = default;

  //! A parameterized constructor
  /*!
    Creates the PolicyRelationship using a list(vector) of relationships.
  */
  explicit PolicyRelationship(const std::vector<GO::Relationship> &relationships) { addRelationships(relationships); }

  //! a method to test if a relatinoship is allowed or not
  /*!
    tests if the relationship is allowed. Overridden to fulfill the PolicyRelationshipInterface
  */
  [[nodiscard]] bool isAllowed(GO::Relationship relationship) const { return relationship_set_.contains(relationship); }

  //! a method to add a relationship to the set of relationships allowed
  /*!
    adds a relationship to the set of relationships allowed by setting its mapped value to true
  */
  void addRelationships(const std::vector<GO::Relationship>& relationships) {

    for (auto const& relation : relationships) {

      relationship_set_.insert(relation);

    }

  }

  //! a method to determine if the Policy is valid
  /*!
    Must contain both Relationship::IS_A and Relationship::PART_OF.
    Cannot contain Relationship::REL_ERROR.
  */
  [[nodiscard]] bool validPolicy() const  {
    return isAllowed(GO::Relationship::IS_A)
           and isAllowed(GO::Relationship::PART_OF)
           and not isAllowed(GO::Relationship::REL_ERROR);
  }

private:
  //! a map of relationships to bool
  /*! maps a relationship to a bool. Boost unordered map give constant time find. */
  OntologySetType<GO::Relationship> relationship_set_;

};

} // namespace

#endif
