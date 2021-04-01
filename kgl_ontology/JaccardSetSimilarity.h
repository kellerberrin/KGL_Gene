/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef JACCARD_SET_SIMILARITY
#define JACCARD_SET_SIMILARITY

#include <SetUtilities.h>
#include <TermSetSimilarityInterface.h>

#include <boost/unordered_set.hpp>

/*! \class JaccardSetSimilarity
	\brief A class to calculate jaccard similarity between 2 sets.

	This class calculates jaccard set similarity between two sets of terms.

*/
class JaccardSetSimilarity : public TermSetSimilarityInterface {

public:

	//! Constructor
	/*!
		Creates the JaccardSetSimilarity class.
	*/
	JaccardSetSimilarity() = default;
  ~JaccardSetSimilarity() override = default;


	//! A method for calculating term set to term set similarity for GO terms;
	/*!
		This method returns the Jaccard set similarity.
	*/
	[[nodiscard]] double calculateSimilarity(const OntologySetType<std::string> &termsA, const  OntologySetType<std::string> &termsB) const override {

		//return 0 if a set is empty
		if(termsA.empty() or termsB.empty()){

			return 0.0;

		}

		//get iterators
		OntologySetType<std::string> _union = SetUtilities::set_union( termsA, termsB);
		OntologySetType<std::string> _intersection = SetUtilities::set_intersection( termsA, termsB);

		if(_union.empty()){

			return 0.0;

		} else {

			return static_cast<double>(_intersection.size())/static_cast<double>(_union.size());

		}

	}

};
#endif
