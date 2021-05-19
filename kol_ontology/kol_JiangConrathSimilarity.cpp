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

  //if the terms do not exist return 0.0 similarity
  if (not ic_map_ptr_->validateTerms(goTermA, goTermB)) {

    return 0.0;

  }

  double maxIC = ic_map_ptr_->getMaxInformation(goTermA);
  //get the MICA value (zero if no mica term)
  double mica_value = ic_map_ptr_->getMICAinfo(goTermA, goTermB, *graph_ptr_);

  double dist = ic_map_ptr_->getValue(goTermA) + ic_map_ptr_->getValue(goTermB) - (2 * mica_value);

  double gosemsim = 1.0 - std::min(1.0, dist/maxIC);

  double standard = 1.0 - (dist / (2.0 * maxIC));

  return GOSemSimFormula ? gosemsim : standard;

}

