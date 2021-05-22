//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_SetSimilarityJaccard.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term set to term set similarity for GO terms;
/*!
  This method returns the Jaccard set similarity.
  Order of row and column terms important when using Asymmetric similarity caches.
*/
double kol::JaccardSetSimilarity::calculateSimilarity( const OntologySetType<std::string> &row_terms,
                                                       const OntologySetType<std::string> &column_terms) const {

  //return 0 if a set is empty
  if (row_terms.empty() or column_terms.empty()) {

    return 0.0;

  }

  //get iterators
  OntologySetType<std::string> set_union = SetUtilities::setUnion(row_terms, column_terms);
  OntologySetType<std::string> set_intersection = SetUtilities::setIntersection(row_terms, column_terms);

  if (set_union.empty()) {

    return 0.0;

  } else {

    return static_cast<double>(set_intersection.size()) / static_cast<double>(set_union.size());

  }

}
