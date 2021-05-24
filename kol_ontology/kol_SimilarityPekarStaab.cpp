//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_SimilarityPekarStaab.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term-to-term similarity for GO terms using Pekar Staab similarity
/*!
  This method returns the PekarStaab similarity.
*/
double kol::SimilarityPekarStaab::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {
  //if the terms do not exit return 0.0 similarity
  if (not depth_map_ptr_->hasTerm(goTermA) or not depth_map_ptr_->hasTerm(goTermB)) {

    return 0.0;

  }

  //if not from same ontology, return 0;
  if (graph_ptr_->getTermOntology(goTermA) != graph_ptr_->getTermOntology(goTermB)) {

    return 0.0;

  }

  //Create 2 sets self + ancestors.
  OntologySetType<std::string> ancestorsA = graph_ptr_->getSelfAncestorTerms(goTermA);
  OntologySetType<std::string> ancestorsB = graph_ptr_->getSelfAncestorTerms(goTermB);

  //if either set is empty, return 0
  if (ancestorsA.empty() || ancestorsB.empty()) {

    return 0.0;

  }

  std::string lca = depth_map_ptr_->getLCA(ancestorsA, ancestorsB);
  size_t lcaDepth = depth_map_ptr_->getValue(lca);
  size_t denom = (depth_map_ptr_->getValue(goTermA) - lcaDepth) + (depth_map_ptr_->getValue(goTermB) - lcaDepth) + lcaDepth;

  if (denom == 0) {

    return 0.0;

  } else {

    return static_cast<double>(lcaDepth) / static_cast<double>(denom);

  }

}

