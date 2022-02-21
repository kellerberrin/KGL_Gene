//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_SetSimilarityAllPairsMax.h"
#include "contrib/kol_Accumulators.h"


namespace kol = kellerberrin::ontology;

//! A method for calculating term set to term set similarity for GO terms;
/*!
  Order of row and column terms important when using Asymmetric similarity caches.
  This method returns the all pairs max similarity.
*/
double kol::SetSimilarityAllPairsMax::calculateSimilarity(const OntologySetType<std::string> &row_terms,
                                                          const OntologySetType<std::string> &column_terms) const {

  //return 0 if a set is empty
  if (row_terms.empty() or column_terms.empty()) {

    return 0.0;

  }

  //get mean accumulator
  Accumulators::MaxAccumulator simMax;

  //iterate A set
  for (auto const &row : row_terms) {
    //iterate B terms
    for (auto const &column : column_terms) {

      //get the term from B set
      double sim = similarity_ptr_->calculateTermSimilarity(row, column);

      //add to accumulator
      simMax(sim);

    }

  }
  //return the mean from the accumulator
  return Accumulators::extractMax(simMax);

}
