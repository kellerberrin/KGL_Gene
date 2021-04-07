/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef PESQUITA_SIMGIC_SET_SIMILARITY
#define PESQUITA_SIMGIC_SET_SIMILARITY

#include <TermSetSimilarityInterface.h>
#include <TermInformationContentMap.h>
#include <SetUtilities.h>
#include <GoGraph.h>


/*! \class PesquitaSimGICSetSimilarity
\brief A class to calculate Pesquita's SimGIC similarity between sets of go terms.

Pesquita, C., Faria, D., Bastos, H., Falcao, A., & Couto, F. (2007, July).
Evaluating GO-based semantic similarity measures. In Proc. 10th Annual
Bio-Ontologies Meeting (Vol. 37, No. 40, p. 38).

*/
class PesquitaSimGICSetSimilarity : public TermSetSimilarityInterface{

public:
	//! Constructor
	/*!
	Creates the PesquitaSimGICSetSimilarity class assigning the GoGraph private member.
	*/
	PesquitaSimGICSetSimilarity(const std::shared_ptr<const GoGraph>& graph, const std::shared_ptr<const TermInformationContentMap>& icMap)
	: _graph(graph), _icMap(icMap) {}
  ~PesquitaSimGICSetSimilarity() override = default;

	//! A method for calculating term set to term set similarity for GO terms;
	/*!
	This method returns the best match average similarity.
	*/
	[[nodiscard]] double calculateSimilarity(const OntologySetType<std::string> &termsA, const OntologySetType<std::string> &termsB) const override {

		// Get the induced set of terms for each set
		OntologySetType<std::string> inducedTermSetA = _graph->getExtendedTermSet(termsA);
    OntologySetType<std::string> inducedTermSetB = _graph->getExtendedTermSet(termsB);
		// Calculate union and intersection
    OntologySetType<std::string> union_set = SetUtilities::set_union(inducedTermSetA, inducedTermSetB);
    OntologySetType<std::string> intersection_set = SetUtilities::set_intersection(inducedTermSetA, inducedTermSetB);

    // Calculate sum for the union set
    double union_sum{0.0};
		for (auto const& term : union_set) {

			union_sum += _icMap->getValue(term);

		}

		// Calculate sum for the intersection set
    double intersection_sum{0.0};
		for (auto const& term : intersection_set) {

			intersection_sum += _icMap->getValue(term);

		}

		//if the union is 0, return 0. No division by 0.
		if (union_sum == 0.0){

			return 0.0;

		} else {

			return intersection_sum / union_sum;

		}

	}

private:

	//! Pointer to the GoGraph object
	/*!
	A reference to GO graph to be used.
	*/
	std::shared_ptr<const GoGraph> _graph;

	//! The information content map
	/*!
	An information content map
	*/
  std::shared_ptr<const TermInformationContentMap> _icMap;

};
#endif
