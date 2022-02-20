//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_SetSimilarityGentlemanSimUI.h"
#include "kol_GoGraphImpl.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term set to term set similarity for GO terms;
/*!
  This method returns the best match average similarity.
  Order of row and column terms important when using Asymmetric similarity caches.
*/


double kol::GentlemanSimUISetSimilarity::calculateSimilarity( const OntologySetType<std::string> &row_terms,
                                                              const OntologySetType<std::string> &column_terms) const {
  // Get the induced set of terms for each set

  OntologySetType<std::string> induced_row_set = graph_ptr_->getGoGraphImpl().getExtendedTermSet(row_terms);
  OntologySetType<std::string> induced_column_set = graph_ptr_->getGoGraphImpl().getExtendedTermSet(column_terms);
  // Calculate union and intersection
  OntologySetType<std::string> union_set = SetUtilities::setUnion(induced_row_set, induced_column_set);
  OntologySetType<std::string> intersection_set = SetUtilities::setIntersection(induced_row_set, induced_column_set);

  //if the union is 0, return 0. No division by 0.
  if (union_set.size() == 0) {

    return 0.0;

  } else {

    return static_cast<double>(intersection_set.size()) / static_cast<double>(union_set.size());

  }

}
