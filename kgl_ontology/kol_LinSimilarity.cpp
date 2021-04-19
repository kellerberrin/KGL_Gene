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

  //if the terms do not exit return 0.0 similarity
  if (not ic_map_ptr_->hasTerm(goTermA) || not ic_map_ptr_->hasTerm(goTermB)) {

    return 0.0;

  }
  //if not from same ontology, return 0;
  if (graph_ptr_->getTermOntology(goTermA) != graph_ptr_->getTermOntology(goTermB)) {

    return 0.0;

  }

  //create 2 sets
  OntologySetType<std::string> ancestorsA = graph_ptr_->getAncestorTerms(goTermA);
  ancestorsA.insert(goTermA);
  OntologySetType<std::string> ancestorsB = graph_ptr_->getAncestorTerms(goTermB);
  ancestorsB.insert(goTermB);

  //if either set is empty, return 0
  if (ancestorsA.empty() or ancestorsB.empty()) {

    return 0.0;

  }

  double mica_info = ic_map_ptr_->getMICAinfo(ancestorsA, ancestorsB);
  //return the normalized information content similarity of Lin
  return (2.0 * mica_info) / (ic_map_ptr_->getValue(goTermA) + ic_map_ptr_->getValue(goTermB));


}
