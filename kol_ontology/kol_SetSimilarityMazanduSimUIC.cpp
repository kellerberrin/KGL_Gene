//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_SetSimilarityMazanduSimUIC.h"
#include "kol_SetUtilities.h"
#include "kol_Accumulators.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term set to term set similarity for GO terms;
/*!
This method returns the best match average similarity.
*/

double kol::SetSimilarityMazanduSimUIC::calculateSimilarity(const OntologySetType<std::string> &row_terms,
                                                            const OntologySetType<std::string> &column_terms) const {

  // Get the induced set of terms for each set
  OntologySetType<std::string> induced_row_terms = graph_ptr_->getExtendedTermSet(row_terms);
  OntologySetType<std::string> induced_column_terms = graph_ptr_->getExtendedTermSet(column_terms);
  // Calculate union and intersection
  OntologySetType<std::string> intersection_set = SetUtilities::setIntersection(induced_row_terms, induced_column_terms);

  double intersection_sum = calcICSum(intersection_set);
  double setA_sum = calcICSum(induced_row_terms);
  double setB_sum = calcICSum(induced_column_terms);


  //if the union is 0, return 0. No division by 0.
  if (setA_sum + setB_sum == 0.0) {

    return 0.0;

  } else {

    if (setA_sum > setB_sum) {

      return intersection_sum / setA_sum;

    } else {

      return intersection_sum / setB_sum;

    }

  }

}

//! A method for calculating the sum of information content of the terms in a set.
/*!
  This method calculates the sum of information content of the terms in a set.
*/

double kol::SetSimilarityMazanduSimUIC::calcICSum(const OntologySetType<std::string> &terms) const {

  double sum{0.0};
  for (auto const &term : terms) {

    sum += ic_map_ptr_->termInformation(term);

  }

  return sum;

}
