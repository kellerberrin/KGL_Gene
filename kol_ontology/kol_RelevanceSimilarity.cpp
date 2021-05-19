//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_RelevanceSimilarity.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term-to-term similarity for GO terms using Relevance similarity
/*!
  This method returns the Relevance similarity.
*/
double kol::RelevanceSimilarity::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {
  //if the terms do not exit return 0.0 similarity
  if (not ic_map_ptr_->validateTerms(goTermA, goTermB)) {

    return 0.0;

  }

  //get the MICA value (zero if no mica term)
  double mica_value = ic_map_ptr_->getMICAinfo(goTermA, goTermB, *graph_ptr_);
  double complement_prob_mica = 1.0 - std::exp(-1.0 * mica_value);
  double denom = (ic_map_ptr_->getValue(goTermA) + ic_map_ptr_->getValue(goTermB));

  if (denom == 0.0 or mica_value == 0.0) {

    return 0.0;

  }

  //return the normalized information content similarity of Relevance
  return ((2.0 * mica_value) / denom) * complement_prob_mica;

}
