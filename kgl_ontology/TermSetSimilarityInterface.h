/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef TERM_SET_SIMILARITY_INTERFACE
#define TERM_SET_SIMILARITY_INTERFACE

#include <string>
#include <boost/unordered_set.hpp>

/*! \class TermSetSimilarityInterface
	\brief An interface class for comparing semantic similarity of sets of GO terms

	This class defines the interface for comparing term set to term set similarity.
	 This is the most useful case for comparing genes with multiple annotations to each other.

*/
class TermSetSimilarityInterface {

public:

  TermSetSimilarityInterface() = default;
  virtual ~TermSetSimilarityInterface() = default;

	//! A pure virtual method for calculating term set to term set similarity for sets of GO terms
	/*!
		This pure virtual method requires similarity measure to implement the basic interface
		  that returns a similarity value for two sets of go terms.
	*/
	[[nodiscard]] virtual double calculateSimilarity(const OntologySetType<std::string> &termsA, const OntologySetType<std::string> &termsB) const = 0;

};
#endif
