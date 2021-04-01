/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef ALL_PAIRS_AVERAGE_SIMILARITY
#define ALL_PAIRS_AVERAGE_SIMILARITY

#include <TermSetSimilarityInterface.h>
#include <TermSimilarityInterface.h>
#include <Accumulators.h>

/*! \class AllPairsAverageSetSimilarity
	\brief A class to calculate the average similarity between all pairs of go terms for 2 sets.

	This class defines the all pairs average similarity getween two sets of terms.
	Put forth by Lord et al.

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, "Investigating semantic similarity
	measures across the Gene Ontology: the relationship between sequence and annotation,"
	Bioinformatics, vol. 19, pp. 1275-83, Jul 1 2003.

*/
class AllPairsAverageSetSimilarity : public TermSetSimilarityInterface{
public:
	//! Constructor
	/*!
		Creates the AllPairsAverageSetSimilarity class assigning the similarity measure private member.
	*/
	AllPairsAverageSetSimilarity(std::shared_ptr<const TermSimilarityInterface> simMeasure) : _similarity(std::move(simMeasure)) {}
  ~AllPairsAverageSetSimilarity() override = default;

	//! A method for calculating term set to term set similarity for GO terms;
	/*!
		This method returns the Relevance similarity.
	*/
	[[nodiscard]] double calculateSimilarity(const OntologySetType<std::string> &termsA, const OntologySetType<std::string> &termsB) const override {
		//return 0 if a set is empty
		if(termsA.empty() or termsB.empty()) {

			return 0.0;

		}

		//get mean accumulator
		Accumulators::MeanAccumulator simMean;

		//get iterators
		//iterate A set
		for(auto const& aTerm : termsA) {
			//iterate B terms
			for(auto const& bTerm : termsB){
				//get the term from B set
				double sim = _similarity->calculateNormalizedTermSimilarity(aTerm,bTerm);
				//std::cout << aTerm << " " << bTerm << " " << sim << std::endl;
				//add to accumulator
				simMean(sim);

			}

		}
		//return the mean from the accumulator
		return Accumulators::extractMean(simMean);

	}

private:
	//! Pointer to the actual similarity measure
	/*!
		This object will actually calculate the similarity between pairs of terms.
	*/
	std::shared_ptr<const TermSimilarityInterface> _similarity;

};
#endif
