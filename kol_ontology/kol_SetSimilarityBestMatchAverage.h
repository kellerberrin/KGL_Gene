/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_BEST_MATCH_AVERAGE_SIMILARITY
#define KGL_BEST_MATCH_AVERAGE_SIMILARITY

#include "kol_SetSimilarityInterface.h"
#include "kol_SimilarityInterface.h"
#include "kol_Accumulators.h"


namespace kellerberrin::ontology {

/*! \class SetSimilarityBestMatchAverage
	\brief A class to calculate the best match average similarity between go terms for 2 sets.

	This class defines the best match average similarity getween two sets of terms.
	Used by Couto et al.

	F. M. Couto, M. J. Silva, and P. M. Coutinho, "Measuring semantic similarity
	between Gene Ontology terms," Data & Knowledge Engineering, vol. 61, 
	pp. 137-152, Apr 2007.
*/
class SetSimilarityBestMatchAverage : public SetSimilarityInterface {

public:
  //! Constructor
  /*!
    Creates the SetSimilarityBestMatchAverage class assigning the similarity measure private memeber.
  */
  explicit SetSimilarityBestMatchAverage(const std::shared_ptr<const SimilarityInterface> &similarity_ptr)
      : similarity_ptr_(similarity_ptr) {}

  ~SetSimilarityBestMatchAverage() override = default;

  //! A method for calculating term set to term set similarity for GO terms;
  /*!
    This method returns the best match average similarity.
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
