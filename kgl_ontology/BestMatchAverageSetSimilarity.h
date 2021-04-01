/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef BEST_MATCH_AVERAGE_SIMILARITY
#define BEST_MATCH_AVERAGE_SIMILARITY

#include <TermSetSimilarityInterface.h>
#include <TermSimilarityInterface.h>
#include <Accumulators.h>

/*! \class BestMatchAverageSetSimilarity
	\brief A class to calculate the best match average similarity between go terms for 2 sets.

	This class defines the best match average similarity getween two sets of terms.
	Used by Couto et al.

	F. M. Couto, M. J. Silva, and P. M. Coutinho, "Measuring semantic similarity
	between Gene Ontology terms," Data & Knowledge Engineering, vol. 61, 
	pp. 137-152, Apr 2007.
*/
class BestMatchAverageSetSimilarity : public TermSetSimilarityInterface{

public:
	//! Constructor
	/*!
		Creates the BestMatchAverageSetSimilarity class assigning the similarity measure private memeber.
	*/
	BestMatchAverageSetSimilarity(std::shared_ptr<const TermSimilarityInterface> simMeasure) : _similarity(std::move(simMeasure)) {}
  ~BestMatchAverageSetSimilarity() override = default;

	//! A method for calculating term set to term set similarity for GO terms;
	/*!
		This method returns the best match average similarity.
	*/
	[[nodiscard]] double calculateSimilarity(const OntologySetType<std::string> &termsA, const OntologySetType<std::string> &termsB) const override {

		//return 0 if a set is empty
		if (termsA.empty() or termsB.empty()) {

			return 0.0;

		}

		//get mean accumulator
		Accumulators::MeanAccumulator simMean;

		//average for best matches
		Accumulators::MeanAccumulator meanA;

		//Have to calculate the best match for term in A to terms in B
		//get iterators
		//iterate A set
		for(auto const& aTerm : termsA) {
			//get best match value
			Accumulators::MaxAccumulator maxForTermA;
			//iterate B terms
			for(auto const bTerm : termsB) {
				//get the term from B set
				double sim = _similarity->calculateNormalizedTermSimilarity( aTerm, bTerm);
				//std::cout << aTerm << " " << bTerm << " " << sim << std::endl;

				//add to accumulator
				maxForTermA(sim);

			}
			//add to accumulator
			meanA(Accumulators::extractMax(maxForTermA));

		}

		//std::cout << "--" << std::endl;

		//Have to calculate the best match for term in B to terms in A
		// then take the average of both so the relationship is symmetric
		//average for best matches
		Accumulators::MeanAccumulator meanB;
		//iterate A set
		for(auto const& bTerm : termsB){
			//get best match value
			Accumulators::MaxAccumulator maxForTermB;
			//iterate B terms
			for(auto const& aTerm : termsA) {
				//get the term from B set
				double sim = _similarity->calculateNormalizedTermSimilarity( aTerm, bTerm);
				//std::cout << aTerm << " " << bTerm << " " << sim << std::endl;

				//add to accumulator
				maxForTermB(sim);

			}
			//add to accumulator
			meanB(Accumulators::extractMax(maxForTermB));

		}

    //Return the average of the 2 means from our accumulator
		double mean_average = (Accumulators::extractMean(meanA) + Accumulators::extractMean(meanB)) / 2.0;

		return mean_average;

	}

private:
	//! Pointer to the actual similarity measure
	/*!
		This object will actually calculate the similarity between pairs of terms.
	*/
	std::shared_ptr<const TermSimilarityInterface> _similarity;

};
#endif
