/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_ALL_PAIRS_AVERAGE_SIMILARITY
#define KOL_ALL_PAIRS_AVERAGE_SIMILARITY

#include "kol_SetSimilarityInterface.h"
#include "kol_SimilarityInterface.h"
#include <memory>

namespace kellerberrin::ontology {

/*! \class SetSimilarityAllPairsAverage
	\brief A class to calculate the average similarity between all pairs of go terms for 2 sets.

	This class defines the all pairs average similarity getween two sets of terms.
	Put forth by Lord et al.

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, "Investigating semantic similarity
	measures across the Gene Ontology: the relationship between sequence and annotation,"
	Bioinformatics, vol. 19, pp. 1275-83, Jul 1 2003.

*/
class SetSimilarityAllPairsAverage : public SetSimilarityInterface {
public:
  //! Constructor
  /*!
    Creates the SetSimilarityAllPairsAverage class assigning the similarity measure private member.
  */
  SetSimilarityAllPairsAverage(const std::shared_ptr<const SimilarityInterface> &similarity_ptr)
      : similarity_ptr_(similarity_ptr) {}

  ~SetSimilarityAllPairsAverage() override = default;

  //! A method for calculating term set to term set similarity for GO terms;
  /*!
    This method returns the Relevance similarity.
    Order of row and column terms important when using Asymmetric similarity caches.
  */
  [[nodiscard]] double calculateSimilarity( const OntologySetType<std::string> &row_terms,
                                            const OntologySetType<std::string> &column_terms) const override;

private:
  //! Pointer to the actual similarity measure
  /*!
    This object will actually calculate the similarity between pairs of terms.
  */
  std::shared_ptr<const SimilarityInterface> similarity_ptr_;

};

} // namespace

#endif
