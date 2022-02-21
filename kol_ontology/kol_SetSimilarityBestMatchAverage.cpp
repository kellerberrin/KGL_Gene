//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_SetSimilarityBestMatchAverage.h"
#include "contrib/kol_Accumulators.h"


namespace kol = kellerberrin::ontology;


//! A method for calculating term set to term set similarity for GO terms;
/*!
  This method returns the best match average similarity.
  Order of row and column terms important when using Asymmetric similarity caches.
*/
double kol::SetSimilarityBestMatchAverage::calculateSimilarity(const OntologySetType<std::string> &row_terms,
                                                               const OntologySetType<std::string> &column_terms) const {

  //return 0 if either set is empty
  if (row_terms.empty() or column_terms.empty()) {

    return 0.0;

  }

  //get mean accumulator
  Accumulators::MeanAccumulator simMean;

  //average for best matches
  Accumulators::MeanAccumulator mean_max_row;

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
    mean_max_row(Accumulators::extractMax(max_row_term));

  }


  //Have to calculate the best match for term in columns to terms in rows
  // then take the average of both so the relationship is symmetric
  //average for best matches
  Accumulators::MeanAccumulator mean_max_column;
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
    mean_max_column(Accumulators::extractMax(max_column_term));

  }

  //Return the average of the 2 means from our accumulator
  double mean_average = (Accumulators::extractMean(mean_max_row) + Accumulators::extractMean(mean_max_column)) / 2.0;

  return mean_average;

}




