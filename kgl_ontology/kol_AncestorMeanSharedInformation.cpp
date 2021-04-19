//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_AncestorMeanSharedInformation.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating the shared infromation between two concepts.
/*!
  This method returns the shared information between two concepts.
*/
double kol::AncestorMeanSharedInformation::sharedInformation(const std::string &termA, const std::string &termB) const  {
  // return 0 for any terms not in the datbase
  if (not ic_map_ptr_->hasTerm(termA) or not ic_map_ptr_->hasTerm(termB)) {

    return 0.0;

  }
  // return 0 for terms in different ontologies
  if (not isSameOntology(termA, termB)) {

    return 0.0;

  }

  Accumulators::MeanAccumulator meanIC;

  OntologySetType<std::string> ancestorsA = graph_ptr_->getSelfAncestorTerms(termA);
  OntologySetType<std::string> ancestorsB = graph_ptr_->getSelfAncestorTerms(termB);

  OntologySetType<std::string> sharedAncestors = SetUtilities::setIntersection(ancestorsA, ancestorsB);

  for (auto const &term : sharedAncestors) {

    meanIC(ic_map_ptr_->getValue(term));

  }

  return Accumulators::extractMean(meanIC);

}

//! An interface method for returning the shared information of a single terms,or information content
/*!
  This method privdes a mechanism for returing a term's infromation content.
*/
double kol::AncestorMeanSharedInformation::sharedInformation(const std::string &term) const {
  // return 0 for any terms not in the datbase
  if (not ic_map_ptr_->hasTerm(term)) {

    return 0.0;

  }

  return ic_map_ptr_->getValue(term);

}


//! An interface method for returning the maximum information content for a term
/*!
  This method provides the absolute max information content within a corpus for normalization purposes.
*/
double kol::AncestorMeanSharedInformation::maxInformationContent(const std::string &term) const {


  //select the correct ontology normalization factor
  GO::Ontology ontology = graph_ptr_->getTermOntology(term);
  double maxIC;

  switch (ontology) {

    case GO::Ontology::BIOLOGICAL_PROCESS:
      maxIC = ic_map_ptr_->getMinBP();
      break;

    case GO::Ontology::MOLECULAR_FUNCTION:
      maxIC = ic_map_ptr_->getMinMF();
      break;

    case GO::Ontology::CELLULAR_COMPONENT:
      maxIC = ic_map_ptr_->getMinCC();
      break;

    default:
    case GO::Ontology::ONTO_ERROR:
      maxIC = 0.0;
      break;

  }

  if (maxIC <= 0.0) {

    return 0.0;

  }

  return -1.0 * std::log(maxIC);

}
