/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef GENTLEMAN_SIMUI_SET_SIMILARITY
#define GENTLEMAN_SIMUI_SET_SIMILARITY

#include <TermSetSimilarityInterface.h>
#include <SetUtilities.h>
#include <GoGraph.h>


/*! \class GentlemanSimUISetSimilarity
	\brief A class to calculate Gentleman's UI similarity between go terms for 2 sets.

	Gentlman R. Visualizing and Distances Using GO. URL http://www.bioconductor.org/docs/vignettes.html.

*/
class GentlemanSimUISetSimilarity : public TermSetSimilarityInterface{

public:
	//! Constructor
	/*!
		Creates the GentlemanUISimilarity class assigning the GoGraph private memeber.
	*/
	GentlemanSimUISetSimilarity(std::shared_ptr<const GoGraph> graph) : _graph(std::move(graph)) {}
	~GentlemanSimUISetSimilarity() override = default;


	//! A method for calculating term set to term set similarity for GO terms;
	/*!
		This method returns the best match average similarity.
	*/
	[[nodiscard]] double calculateSimilarity(const OntologySetType<std::string> &termsA, const OntologySetType<std::string> &termsB ) const override {
		// Get the induced set of terms for each set

		OntologySetType<std::string> inducedTermSetA = _graph->getExtendedTermSet(termsA);
		OntologySetType<std::string> inducedTermSetB = _graph->getExtendedTermSet(termsB);
		// Calculate union and intersection
		OntologySetType<std::string> union_set = SetUtilities::set_union(inducedTermSetA, inducedTermSetB);
		OntologySetType<std::string> intersection_set = SetUtilities::set_intersection(inducedTermSetA, inducedTermSetB);

		//if the union is 0, return 0. No division by 0.
		if (union_set.size() == 0){

			return 0.0;

		} else {

			return static_cast<double>(intersection_set.size()) / static_cast<double>(union_set.size());

		}

	}

private:

	//! Pointer to the GoGraph object
	/*!
		A reference to GO graph to be used.
	*/
	std::shared_ptr<const GoGraph> _graph;



};
#endif
