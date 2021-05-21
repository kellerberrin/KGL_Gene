//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_AllPairsAverageSetSimilarity.h"

namespace kol = kellerberrin::ontology;


double kol::AllPairsAverageSetSimilarity::calculateSimilarity( const OntologySetType<std::string> &row_terms,
                                                               const OntologySetType<std::string> &column_terms) const {
  //return 0 if a set is empty
  if (row_terms.empty() or column_terms.empty()) {

    return 0.0;

  }

  //get mean accumulator
  Accumulators::MeanAccumulator simMean;

  //get iterators
  //iterate A set
  for (auto const &row : row_terms) {
    //iterate B terms
    for (auto const &column : column_terms) {

      //get the term from B set
      double sim = similarity_ptr_->calculateTermSimilarity(row, column);
      //add to accumulator
      simMean(sim);

    }

  }
  //return the mean from the accumulator
  return Accumulators::extractMean(simMean);

}
