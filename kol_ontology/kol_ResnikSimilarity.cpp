//
// Created by kellerberrin on 19/4/21.
//

#include "kol_ResnikSimilarity.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term-to-term similarity for GO terms using Resnik similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/
double kol::ResnikSimilarity::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {
  //if the terms do not exit return 0.0 similarity
  if (not ic_map_ptr_->validateTerms(goTermA, goTermB)) {

    return 0.0;

  }

  return ic_map_ptr_->getMICAinfo(goTermA, goTermB, *graph_ptr_);

}

//! A method for calculating term-to-term similarity for GO terms using Normalized Resnik similarity
/*!
  This method returns the Resnik similarity divided by the maximum possible similarity
*/
double kol::ResnikSimilarity::calculateNormalizedTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {
  //call base similarity
  double resnik = calculateTermSimilarity(goTermA, goTermB);

  if (resnik <= 0) {

    return 0.0;

  }

  double maxIC = ic_map_ptr_->getMaxInformation(goTermA);
  return resnik / maxIC;

}
