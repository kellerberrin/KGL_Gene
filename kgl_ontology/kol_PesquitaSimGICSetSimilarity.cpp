//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_PesquitaSimGICSetSimilarity.h"

namespace kol = kellerberrin::ontology;


[[nodiscard]] double kol::PesquitaSimGICSetSimilarity::calculateSimilarity(const OntologySetType<std::string> &row_terms,
                                                                           const OntologySetType<std::string> &column_terms) const {

  // Get the induced set of terms for each set
  OntologySetType<std::string> induced_row_terms = graph_ptr_->getExtendedTermSet(row_terms);
  OntologySetType<std::string> induced_column_terms = graph_ptr_->getExtendedTermSet(column_terms);
  // Calculate union and intersection
  OntologySetType<std::string> union_set = SetUtilities::setUnion(induced_row_terms, induced_column_terms);
  OntologySetType<std::string> intersection_set = SetUtilities::setIntersection(induced_row_terms, induced_column_terms);

  // Calculate sum for the union set
  double union_sum{0.0};
  for (auto const &term : union_set) {

    union_sum += ic_map_ptr_->getValue(term);

  }

  // Calculate sum for the intersection set
  double intersection_sum{0.0};
  for (auto const &term : intersection_set) {

    intersection_sum += ic_map_ptr_->getValue(term);

  }

  //if the union is 0, return 0. No division by 0.
  if (union_sum == 0.0) {

    return 0.0;

  } else {

    return intersection_sum / union_sum;

  }

}
