//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_SetSimilarityMazanduSimDIC.h"
#include <kol_SetUtilities.h>

namespace kol = kellerberrin::ontology;


double kol::SetSimilarityMazanduSimDIC::calculateSimilarity(const OntologySetType<std::string> &row_terms,
                                                            const OntologySetType<std::string> &column_terms) const {

  // Get the induced set of terms for each set
  OntologySetType<std::string> induced_row_terms = graph_->getExtendedTermSet(row_terms);
  OntologySetType<std::string> induced_column_terms = graph_->getExtendedTermSet(column_terms);
  // Calculate union and intersection
  OntologySetType<std::string> intersection_set = SetUtilities::setIntersection(induced_row_terms, induced_column_terms);

  double intersection_sum = calcICSum(intersection_set);
  double setA_sum = calcICSum(induced_row_terms);
  double setB_sum = calcICSum(induced_column_terms);

  //if the union is 0, return 0. No division by 0.
  if (setA_sum + setB_sum == 0.0) {

    return 0.0;

  } else {

    return 2 * intersection_sum / (setA_sum + setB_sum);

  }

}


double kol::SetSimilarityMazanduSimDIC::calcICSum(const OntologySetType<std::string> &terms) const {

  double sum{0.0};
  for (auto const &term : terms) {

    sum += info_content_map_->termInformation(term);

  }

  return sum;

}
