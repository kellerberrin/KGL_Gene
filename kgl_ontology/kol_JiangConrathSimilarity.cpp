//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_JiangConrathSimilarity.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term-to-term similarity for GO terms using JiangConrath similarity
/*!
  This method returns the JiangConrath similarity or the information content of the most informative common ancestor.
*/
double kol::JiangConrathSimilarity::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {
  //if the terms do not exit return 0.0 similarity
  if (not ic_map_ptr_->hasTerm(goTermA) or not ic_map_ptr_->hasTerm(goTermB)) {

    return 0.0;

  }

  //if not from same ontology, return 0;
  if (graph_ptr_->getTermOntology(goTermA) != graph_ptr_->getTermOntology(goTermB)) {

    return 0;

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

  double maxIC;
  //select the correct ontology normalization factor
  GO::Ontology ontoType = graph_ptr_->getTermOntology(goTermA);
  if (ontoType == GO::Ontology::BIOLOGICAL_PROCESS) {

    maxIC = -std::log(ic_map_ptr_->getMinBP());

  } else if (ontoType == GO::Ontology::MOLECULAR_FUNCTION) {

    maxIC = -std::log(ic_map_ptr_->getMinMF());

  } else {

    maxIC = -std::log(ic_map_ptr_->getMinCC());

  }

  //get the MICA value (zero if no mica term)
  double mica_value = ic_map_ptr_->getMICAinfo(ancestorsA, ancestorsB);

  double dist = ic_map_ptr_->getValue(goTermA) + ic_map_ptr_->getValue(goTermB) - (2 * mica_value);

  return 1.0 - (dist / (2.0 * maxIC));

}

