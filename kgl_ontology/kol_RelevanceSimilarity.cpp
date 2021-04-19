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
  if (not ic_map_ptr_->hasTerm(goTermA) or not ic_map_ptr_->hasTerm(goTermB)) {

    return 0.0;

  }

  //if not from same ontology, return 0;
  if (graph_ptr_->getTermOntology(goTermA) != graph_ptr_->getTermOntology(goTermB)) {

    return 0.0;

  }

  //create 2 sets
  OntologySetType<std::string> ancestorsA = graph_ptr_->getSelfAncestorTerms(goTermA);
  OntologySetType<std::string> ancestorsB = graph_ptr_->getSelfAncestorTerms(goTermB);

  //if either set is empty, return 0
  if (ancestorsA.empty() or ancestorsB.empty()) {

    return 0.0;

  }

  //get the MICA
  double mica_info = ic_map_ptr_->getMICAinfo(ancestorsA, ancestorsB);
  double complement_prob_mica = 1.0 - std::exp(-1.0 * mica_info);
  double denom = (ic_map_ptr_->getValue(goTermA) + ic_map_ptr_->getValue(goTermB));

  if (denom == 0.0 or mica_info == 0.0) {

    return 0.0;

  }

  //return the normalized information content similarity of Relevance
  return ((2.0 * mica_info) / denom) * complement_prob_mica;

}
