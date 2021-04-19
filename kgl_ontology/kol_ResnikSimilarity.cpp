//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_ResnikSimilarity.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term-to-term similarity for GO terms using Resnik similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/
double kol::ResnikSimilarity::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {
  //if the terms do not exit return 0.0 similarity
  if (not ic_map_ptr_->hasTerm(goTermA) or not ic_map_ptr_->hasTerm(goTermB)) {

    return 0.0;

  }

  //if not from same ontology, return 0;
  if (graph_ptr_->getTermOntology(goTermA) != graph_ptr_->getTermOntology(goTermB)) {

    return 0.0;

  }

  //Create 2 sets of term + ancestors
  OntologySetType<std::string> ancestorsA = graph_ptr_->getSelfAncestorTerms(goTermA);
  OntologySetType<std::string> ancestorsB = graph_ptr_->getSelfAncestorTerms(goTermB);

  return ic_map_ptr_->getMICAinfo(ancestorsA, ancestorsB);

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

  //select the correct ontology normalization factor
  GO::Ontology ontology = graph_ptr_->getTermOntology(goTermA);
  double ontology_info;
  switch (ontology) {

    case GO::Ontology::BIOLOGICAL_PROCESS:
      ontology_info = ic_map_ptr_->getMinBP();
      break;

    case GO::Ontology::MOLECULAR_FUNCTION:
      ontology_info = ic_map_ptr_->getMinMF();
      break;

    case GO::Ontology::CELLULAR_COMPONENT:
      ontology_info = ic_map_ptr_->getMinCC();
      break;

    default:
    case GO::Ontology::ONTO_ERROR:
      ontology_info = 0.0;
      break;
  }

  if (ontology_info <= 0.0) {

    return 0.0;

  }

  double maxIC = -1.0 * std::log(ontology_info);
  return resnik / maxIC;

}
