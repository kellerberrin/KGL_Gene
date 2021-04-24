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
  AllPairsMaxSetSimilarity(const std::shared_ptr<const TermSimilarityInterface> &similarity_ptr)
      : similarity_ptr_(similarity_ptr) {}

  ~AllPairsMaxSetSimilarity() override = default;


  //! A method for calculating term set to term set similarity for GO terms;
  /*!
    Order of row and column terms important when using Asymmetric similarity caches.
    This method returns the all pairs max similarity.
  */
  [[nodiscard]] double calculateSimilarity( const OntologySetType<std::string> &row_terms,
                                            const OntologySetType<std::string> &column_terms) const override;

private:
  //! Pointer to the actual similarity measure
  /*!
    This object will actually calculate the similarity between pairs of terms.
  */
  std::shared_ptr<const TermSimilarityInterface> similarity_ptr_;

};

} // namespace

#endif
