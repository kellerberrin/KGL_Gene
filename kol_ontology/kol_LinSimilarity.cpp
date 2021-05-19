//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_LinSimilarity.h"

namespace kol = kellerberrin::ontology;


double kol::LinSimilarity::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {

  if (goTermA == goTermB) {

    return 1.0;

  }

  //If the terms do not exit of have different ontology return 0.0
  if (not ic_map_ptr_->validateTerms(goTermA, goTermB)) {

    return 0.0;

  }

  //get the MICA value (zero if no mica term)
  double mica_value = ic_map_ptr_->getMICAinfo(goTermA, goTermB, *graph_ptr_);
  //return the normalized information content similarity of Lin
  return (2.0 * mica_value) / (ic_map_ptr_->getValue(goTermA) + ic_map_ptr_->getValue(goTermB));

}
