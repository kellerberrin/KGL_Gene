/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef MAZANDU_SIMUIC_SET_SIMILARITY
#define MAZANDU_SIMUIC_SET_SIMILARITY

#include <TermSetSimilarityInterface.h>
#include <TermInformationContentMap.h>
#include <SetUtilities.h>
#include <Accumulators.h>
#include <GoGraph.h>


/*! \class MazanduSimUICSetSimilarity
\brief A class to calculate Mazandu and Mulder's SimUIC similarity between 2 sets of go terms.

A separate measure from their SimDIC.

Mazandu, G. K., & Mulder, N. J. (2014). Information content-based Gene Ontology functional
similarity measures: which one to use for a given biological data type?. PloS one, 9(12), e113859

*/
class MazanduSimUICSetSimilarity : public TermSetSimilarityInterface{

public:
	//! Constructor
	/*!
	Creates the MazanduSimUICSetSimilarity class assigning the GoGraph private memeber.
	*/
	MazanduSimUICSetSimilarity(std::shared_ptr<const GoGraph> graph, std::shared_ptr<const TermInformationContentMap> icMap)
	: _graph(std::move(graph)), _icMap(std::move(icMap)) {}
  ~MazanduSimUICSetSimilarity() override = default;


	//! A method for calculating term set to term set similarity for GO terms;
	/*!
	This method returns the best match average similarity.
	*/
	[[nodiscard]] double calculateSimilarity(const OntologySetType<std::string> &termsA, const OntologySetType<std::string> &termsB) const override {

		// Get the induced set of terms for each set
		OntologySetType<std::string> inducedTermSetA = _graph->getExtendedTermSet(termsA);
		OntologySetType<std::string> inducedTermSetB = _graph->getExtendedTermSet(termsB);
		// Calculate union and intersection
    OntologySetType<std::string> intersection_set = SetUtilities::set_intersection(inducedTermSetA, inducedTermSetB);

		double intersection_sum = calcICSum(intersection_set);
		double setA_sum = calcICSum(inducedTermSetA);
		double setB_sum = calcICSum(inducedTermSetB);


		//if the union is 0, return 0. No division by 0.
		if (setA_sum + setB_sum == 0.0) {

			return 0.0;

		} else {

			if (setA_sum > setB_sum){

				return intersection_sum / setA_sum;

			} else{

				return intersection_sum / setB_sum;

			}
			
		}

	}

private:

	//! Pointer to the GoGraph object.
	/*!
	A reference to GO graph to be used.
	*/
	std::shared_ptr<const GoGraph> _graph;

	//! The information content map.
	/*!
	An information content map.
	*/
	std::shared_ptr<const TermInformationContentMap> _icMap;


	//! A method for calculating the sum of information content of the terms in a set.
	/*!
		This method calculates the sum of information content of the terms in a set.
	*/
	[[nodiscard]] double calcICSum(const OntologySetType<std::string> &terms) const {

		double sum{0.0};
		for (auto const& term : terms) {

			sum += _icMap->getValue(term);

		}

		return sum;

	}

};
#endif
