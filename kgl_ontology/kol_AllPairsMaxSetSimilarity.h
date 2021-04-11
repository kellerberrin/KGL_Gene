/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ALL_PAIRS_MAX_SIMILARITY
#define KGL_ALL_PAIRS_MAX_SIMILARITY

#include "kol_TermSetSimilarityInterface.h"
#include "kol_TermSimilarityInterface.h"
#include "kol_Accumulators.h"

namespace kellerberrin::ontology {


/*! \class AllPairsMaxSetSimilarity
	\brief A class to calculate the max similarity between all pairs of go terms for 2 sets.

	This class defines the all pairs max similarity between two sets of terms.
	Used by Sevilla et al.
	
	J. L. Sevilla, V. Segura, A. Podhorski, E. Guruceaga, J. M. Mato, 
	 L. A. Martinez-Cruz, et al., "Correlation between gene expression and
	 GO semantic similarity," IEEE/ACM Trans Comput Biol Bioinform, 
	 vol. 2, pp. 330-8, Oct-Dec 2005.

*/
class AllPairsMaxSetSimilarity : public TermSetSimilarityInterface {


public:

  //! Constructor
  /*!
    Creates the AllPairsMaxSetSimilarity class assigning the similarity measure private member.
  */
  AllPairsMaxSetSimilarity(const std::shared_ptr<const TermSimilarityInterface> &simMeasure)
      : _similarity(simMeasure) {}

  ~AllPairsMaxSetSimilarity() override = default;


  //! A method for calculating term set to term set similarity for GO terms;
  /*!
    This method returns the all pairs max similarity.
  */
  [[nodiscard]] double calculateSimilarity(const OntologySetType<std::string> &termsA, const OntologySetType<std::string> &termsB) const override {

    //return 0 if a set is empty
    if (termsA.empty() or termsB.empty()) {

      return 0.0;

    }

    //get mean accumulator
    Accumulators::MaxAccumulator simMax;

    //iterate A set
    for (auto const &aTerm : termsA) {
      //iterate B terms
      for (auto const &bTerm : termsB) {
        //get the term from B set
        double sim = _similarity->calculateNormalizedTermSimilarity(aTerm, bTerm);
        //std::cout << aTerm << " " << bTerm << " " << sim << std::endl;
        //add to accumulator
        simMax(sim);

      }

    }
    //return the mean from the accumulator
    return Accumulators::extractMax(simMax);

  }

private:
  //! Pointer to the actual similarity measure
  /*!
    This object will actually calculate the similarity between pairs of terms.
  */
  std::shared_ptr<const TermSimilarityInterface> _similarity;

};

} // namespace

#endif
