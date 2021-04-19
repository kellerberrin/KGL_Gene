/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_JACCARD_SET_SIMILARITY
#define KGL_JACCARD_SET_SIMILARITY

#include "kol_SetUtilities.h"
#include "kol_TermSetSimilarityInterface.h"

namespace kellerberrin::ontology {

/*! \class JaccardSetSimilarity
	\brief A class to calculate jaccard similarity between 2 sets.

	This class calculates jaccard set similarity between two sets of terms.

*/
class JaccardSetSimilarity : public TermSetSimilarityInterface {

public:

  //! Constructor
  /*!
    Creates the JaccardSetSimilarity class.
  */
  JaccardSetSimilarity() = default;

  ~JaccardSetSimilarity() override = default;


  //! A method for calculating term set to term set similarity for GO terms;
  /*!
    This method returns the Jaccard set similarity.
    Order of row and column terms important when using Asymmetric similarity caches.
  */
  [[nodiscard]] double calculateSimilarity( const OntologySetType<std::string> &row_terms,
                                            const OntologySetType<std::string> &column_terms) const override;

};

} // namespace

#endif
