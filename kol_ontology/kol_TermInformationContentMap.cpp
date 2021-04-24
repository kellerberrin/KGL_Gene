//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_TermInformationContentMap.h"

namespace kol = kellerberrin::ontology;


//! A parameterized constructor
/*!
  This constructor takes pointers to GoGraph and AnnotationData objects.
    Only the parameterized construtor is allowed to ensure these objects are
    created with valid parameters.
    This constructor relies on the TermProbabilityMap.
*/
kol::TermInformationContentMap::TermInformationContentMap( const std::shared_ptr<const GoGraph> &graph,
                                                           const std::shared_ptr<const AnnotationData> &annoData)
    : TermProbabilityMap(graph, annoData) {

  for (size_t i = 0; i < probabilities().size(); ++i) {


    if (probabilities().at(i) <= 0.0) {

      probabilities().at(i) = BAD_INFO_VALUE_;

    } else {

      probabilities().at(i) = -1.0 * std::log(probabilities().at(i));

    }


  }//end for, each probability value

}
