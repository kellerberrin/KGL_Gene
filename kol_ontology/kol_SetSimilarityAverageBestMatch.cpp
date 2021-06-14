//
// Created by kellerberrin on 10/6/21.
//

#include "kol_SetSimilarityAverageBestMatch.h"
#include "kol_OntologyTypes.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term set to term set similarity for GO terms;
/*!
  This method returns the best match average similarity.
  Order of row and column terms important when using Asymmetric similarity caches.
*/
double kol::SetSimilarityAverageBestMatch::calculateSimilarity(const OntologySetType<std::string> &row_terms,
                                                               const OntologySetType<std::string> &column_terms) const {

  //return 0 if either set is empty
  if (row_terms.empty() or column_terms.empty()) {

    return 0.0;

  }

  double sum_max{0};

  //Have to calculate the best match for term in rows to terms in columns
  //get iterators
  //iterate A set
  for (auto const &row : row_terms) {
    //get best match value
    Accumulators::MaxAccumulator max_row_term;
    //iterate B terms
    for (auto const& column : column_terms) {
      //get the term from B set
      double sim = similarity_ptr_->calculateTermSimilarity(row, column);

      //add to accumulator
      max_row_term(sim);

    }
    //add to accumulator
    sum_max += Accumulators::extractMax(max_row_term);

  }


  //Have to calculate the best match for term in columns to terms in rows
  // then take the average of both so the relationship is symmetric
  //average for best matches
  //iterate A set
  for (auto const &column : column_terms) {
    //get best match value
    Accumulators::MaxAccumulator max_column_term;
    //iterate B terms
    for (auto const &row : row_terms) {
      //get the term from B set
      double sim = similarity_ptr_->calculateTermSimilarity(row, column);

      //add to accumulator
      max_column_term(sim);

    }
    //add to accumulator
    sum_max += Accumulators::extractMax(max_column_term);

  }

  //Return the average of the 2 means from our accumulator
  double mean_average = sum_max / static_cast<double>(column_terms.size() + row_terms.size());

  return mean_average;

}




