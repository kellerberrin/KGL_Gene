/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_TERM_DEPTH_MAP
#define KOL_TERM_DEPTH_MAP

#include "kol_TermAnnotation.h"

#include "kol_GoGraph.h"


namespace kellerberrin::ontology {

/*! \class InformationDepthMap
	\brief A class to calculate the depth of a GO term in the ontology.

	This class provides a map that returns the depth of a GO term. This method
	 is used in graph and edge based similarity methods to calculate a node's depth
*/

class InformationDepthMap {
public:
  //! A parameterized constructor
  /*!
    This constructor takes pointers to GoGraphImpl and TermAnnotation objects.
      Only the parameterized construtor is allowed to ensure these objects are
      created with valid parameters.
  */
  explicit InformationDepthMap(const GoGraphImpl &graph) { initializeDepthMap(graph); }
  //! A default constructor
  /*!
    This constructor initialized the storage structures. Should not be used.
  */
  InformationDepthMap() = delete;
  //! A default desctructor
  /*!
    This desctructor clears the containters
  */
  ~InformationDepthMap() = default;

  //! Accessor for probablities vector
  /*!
    Get the vector of values
  */
  [[nodiscard]] const std::vector<size_t> &getValues() const { return depths_; }

  //! Function to return all the keys in the map
  /*!
    Returns all valid keys in the map.
  */
  [[nodiscard]] std::vector<std::string> getKeys() const;

  //! Method to test if the id exists in the map
  /*!
    Return true the id is found, false if not
  */
  [[nodiscard]] bool hasTerm(const std::string &testTerm) const { return name_to_index_.find(testTerm) != name_to_index_.end(); }

  //! Mapping function to return the value mapped by key
  /*!
    Get the value mapped by the given key. A specified function for the [] operator
  */
  [[nodiscard]] size_t getValue(const std::string &termId) const;

  //! A method for calculating the least common ancestor
  /*!
    This method searches the sets to determine the deepest common ancestor
  */
  [[nodiscard]] std::string getLCA( const OntologySetType<std::string> &ancestorsA,
                                    const OntologySetType<std::string> &ancestorsB) const;

private:
  //! A private map that returns the index of a term.
  /*!
    This map takes string term ids and returns the index for annotation count access.
  */
  OntologyMapType<std::string, size_t> name_to_index_;

  //! A private list of term depths
  /*!
    This vector of doubles holds the depth for each term
  */
  std::vector<size_t> depths_;

  //! A private method to calculate the depth values on object construction
  /*!
    This method actually calculates the depth values.
  */

  void initializeDepthMap(const GoGraphImpl &graph);

  //! A method for calculating the least common ancestor
  /*!
    This method searches the sets to determine the deepest common ancestor
  */
  [[nodiscard]] std::string getEfficientLCA( const OntologySetType<std::string> &smaller_set,
                                             const OntologySetType<std::string> &larger_set) const;

};

} // namespace

#endif
