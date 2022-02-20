//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_InformationAncestorMean.h"
#include "kol_GoGraphImpl.h"

#include "kol_SetUtilities.h"
#include "kol_Accumulators.h"


#include <boost/accumulators/statistics/max.hpp>


namespace kol = kellerberrin::ontology;


//! A method for calculating the shared infromation between two concepts.
/*!
  This method returns the shared information between two concepts.
*/
double kol::InformationAncestorMean::sharedInformation(const std::string &termA, const std::string &termB) const  {
  // return 0 for any terms not in the datbase
  if (not ic_map_ptr_->validateTerms(termA, termB)) {

    return 0.0;

  }

  Accumulators::MeanAccumulator meanIC;

  OntologySetType<std::string> ancestorsA = graph_ptr_->getGoGraphImpl().getSelfAncestorTerms(termA);
  OntologySetType<std::string> ancestorsB = graph_ptr_->getGoGraphImpl().getSelfAncestorTerms(termB);

  OntologySetType<std::string> sharedAncestors = SetUtilities::setIntersection(ancestorsA, ancestorsB);

  for (auto const &term : sharedAncestors) {

    meanIC(ic_map_ptr_->termInformation(term));

  }

  return Accumulators::extractMean(meanIC);

}

