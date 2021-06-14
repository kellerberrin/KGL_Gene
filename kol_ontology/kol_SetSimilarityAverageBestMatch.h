//
// Created by kellerberrin on 10/6/21.
//

#ifndef KOL_SETSIMILARITYBMA_ALT_H
#define KOL_SETSIMILARITYBMA_ALT_H


#include "kol_SetSimilarityInterface.h"
#include "kol_SimilarityInterface.h"
#include "kol_Accumulators.h"


namespace kellerberrin::ontology {

/*! \class SetSimilarityBestMatchAverage
	\brief A class to calculate the Average Best Match similarity between go terms for 2 sets.

	This class defines the Average Best Match (ABM) similarity between two sets of terms.
  Not to be confused with Best Match Average (BMA).

 "DaGO-Fun: tool for Gene Ontology-based functional analysis using term information content measures"
 Gaston K Mazandu and Nicola J Mulder; BMC Bioinformatics 2013, 14:284.

*/
class SetSimilarityAverageBestMatch : public SetSimilarityInterface {

public:
  //! Constructor
  /*!
    Creates the SetSimilarityBestMatchAverage class assigning the similarity measure private memeber.
  */
  explicit SetSimilarityAverageBestMatch(const std::shared_ptr<const SimilarityInterface> &similarity_ptr)
      : similarity_ptr_(similarity_ptr) {}

  ~SetSimilarityAverageBestMatch() override = default;

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


#endif //KGL_KOL_SETSIMILARITYBMA_ALT_H
