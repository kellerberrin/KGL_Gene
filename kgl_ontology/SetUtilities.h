/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef SET_UTILITIES
#define SET_UTILITIES

#include <OntologyTypes.h>
#include <boost/unordered_set.hpp>
#include <boost/foreach.hpp>

//! The SetUtilities namespace provides useful opperations on generic boost::unordered_set
/*!
	SetUtilities provides useful operations on boost::unordered_set. 
	set_intersection is optimized to traverse the smaller set and test exisitence in the larger set.
*/
class SetUtilities{

public:
  // Ony static definitions.
  SetUtilities() = delete;

	//!set UTILITY function to determine if a set contains an element
	template<typename Type>
	[[nodiscard]] static bool set_contains(const OntologySetType<Type> &mySet, const Type& element)
	{

		return mySet.find(element) != mySet.end();

	}

	//!set UTILITY funciton convert a vector to a set
	template<typename Type>
  [[nodiscard]] static OntologySetType<Type> convert_vector(const std::vector<Type> &myVec)
	{

		OntologySetType<Type> s(myVec.begin(), myVec.end());
		return s;

	}

	//!set OPERATION intersection for boost::unordered_set
	template<typename Type>
  [[nodiscard]] static OntologySetType<Type> set_intersection(const OntologySetType<Type> &setA, const OntologySetType<Type> &setB)
	{
		//iterate the smaller set O(n)
		if(setA.size() <= setB.size()){

			OntologySetType<Type> intersectionSet;
			for (auto const &element : setA ) {
				//find element in O(1)
				if(set_contains(setB,element)){

					intersectionSet.insert(element);

				}

			}

			return intersectionSet;

		} else {

			return set_intersection(setB,setA);

		}

	}

	//!set OPERATION union for boost::unordered_set
	template<typename Type>
  [[nodiscard]] static OntologySetType<Type> set_union(const OntologySetType<Type> &setA, const OntologySetType<Type> &setB)
	{

		if(setA.size() <= setB.size()) {

			OntologySetType<Type> unionSet = setB;
			unionSet.insert(setA.begin(),setA.end());
			return unionSet;

		} else {

			return set_union(setB,setA);

		}

	}

	//!set OPERATION difference A - B
	template<typename Type>
  [[nodiscard]] static OntologySetType<Type> set_difference(const OntologySetType<Type> &setA, const OntologySetType<Type> &setB)
	{

		OntologySetType<Type> differenceSet;
		for (auto const &element : setA) {

			if(!set_contains(setB,element)) {

				differenceSet.insert(element);

			}

		}
		return differenceSet;

	}

};

#endif
