#ifndef KOL_SET_UTILITIES
#define KOL_SET_UTILITIES

#include "kol_OntologyTypes.h"

#include <vector>


namespace kellerberrin::ontology {

//! The SetUtilities namespace provides useful operations on generic boost::unordered_set
/*!
	SetUtilities provides useful operations on boost::unordered_set. 
	setIntersection is optimized to traverse the smaller set and test exisitence in the larger set.
*/
class SetUtilities {

public:
  // Ony static definitions.
  SetUtilities() = delete;

  //!set UTILITY function convert a vector to a set
  template<typename Type>
  [[nodiscard]] static OntologySetType<Type> convertVector(const std::vector<Type> &vector) {

    OntologySetType<Type> set(vector.begin(), vector.end());
    return set;

  }

  //!set UTILITY function convert a set to a vector
  template<typename Type>
  [[nodiscard]] static std::vector<Type> convertSet(const OntologySetType<Type>& set) {

    std::vector<Type> vector(set.begin(), set.end());
    return vector;

  }

  //!set OPERATION intersection f
  template<typename Type>
  [[nodiscard]] static OntologySetType<Type> setIntersection(const OntologySetType<Type> &setA, const OntologySetType<Type> &setB) {
    //iterate the smaller set O(n)
    if (setA.size() <= setB.size()) {

      OntologySetType<Type> intersectionSet;
      for (auto const &element : setA) {
        //find element in O(1)
        if (setB.contains(element)) {

          intersectionSet.insert(element);

        }

      }

      return intersectionSet;

    } else {

      return setIntersection(setB, setA);

    }

  }

  //!set OPERATION union
  template<typename Type>
  [[nodiscard]] static OntologySetType<Type> setUnion(const OntologySetType<Type> &setA, const OntologySetType<Type> &setB) {

    if (setA.size() <= setB.size()) {

      OntologySetType<Type> unionSet = setB;
      unionSet.insert(setA.begin(), setA.end());
      return unionSet;

    } else {

      return setUnion(setB, setA);

    }

  }

  //!set OPERATION difference A - B
  template<typename Type>
  [[nodiscard]] static OntologySetType<Type> setDifference(const OntologySetType<Type> &setA, const OntologySetType<Type> &setB) {

    OntologySetType<Type> differenceSet;
    for (auto const &element : setA) {

      if (not setB.contains(element)) {

        differenceSet.insert(element);

      }

    }
    return differenceSet;

  }

};

} // namespace

#endif
