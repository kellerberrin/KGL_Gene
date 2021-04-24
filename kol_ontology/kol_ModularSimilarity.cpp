//
// Created by kellerberrin on 19/4/21.
//

#include "kol_ModularJiangConrath.h"
#include "kol_ModularLin.h"
#include "kol_ModularResnik.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term-to-term similarity for GO terms using JiangConrath similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/

double kol::ModularJiangConrath::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {

  if (goTermA == goTermB) {

    return 1.0;

  }

  if (not shared_info_ptr_->hasTerm(goTermA) or not shared_info_ptr_->hasTerm(goTermB)) {

    return 0.0;

  }

  if (not shared_info_ptr_->isSameOntology(goTermA, goTermB)) {

    return 0.0;

  }

  double sharedIC = shared_info_ptr_->sharedInformation(goTermA, goTermB);
  double termA_IC = shared_info_ptr_->sharedInformation(goTermA);
  double termB_IC = shared_info_ptr_->sharedInformation(goTermB);
  double maxIC = shared_info_ptr_->maxInformationContent(goTermA);

  double dist = termA_IC + termB_IC - (2.0 * sharedIC);

  return 1.0 - (dist / (2.0 * maxIC));

}


//! A method for calculating term-to-term similarity for GO terms using Lin similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/
double kol::ModularLin::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {

  if (goTermA == goTermB) {

    return 1.0;

  }

  if (not shared_info_ptr_->hasTerm(goTermA) or not shared_info_ptr_->hasTerm(goTermB)) {

    return 0.0;

  }

  if (not shared_info_ptr_->isSameOntology(goTermA, goTermB)) {

    return 0.0;

  }

  double sharedIC = shared_info_ptr_->sharedInformation(goTermA, goTermB);
  double termA_IC = shared_info_ptr_->sharedInformation(goTermA);
  double termB_IC = shared_info_ptr_->sharedInformation(goTermB);

  return (2.0 * sharedIC) / (termA_IC + termB_IC);

}


//! A method for calculating term-to-term similarity for GO terms using Resnik similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/
double kol::ModularResnik::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {

  if (not shared_info_ptr_->hasTerm(goTermA) or not shared_info_ptr_->hasTerm(goTermB)) {

    return 0.0;

  }
  if (not shared_info_ptr_->isSameOntology(goTermA, goTermB)) {

    return 0.0;

  }

  return shared_info_ptr_->sharedInformation(goTermA, goTermB);

}


//! A method for calculating term-to-term similarity for GO terms using Normalized Resnik similarity
/*!
  This method returns the Resnik similarity divided by the maximum possible similarity
*/
double kol::ModularResnik::calculateNormalizedTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {

  if (not shared_info_ptr_->hasTerm(goTermA) or not shared_info_ptr_->hasTerm(goTermB)) {

    return 0.0;

  }
  if (not shared_info_ptr_->isSameOntology(goTermA, goTermB)) {

    return 0.0;

  }

  double sharedInformation = shared_info_ptr_->sharedInformation(goTermA, goTermB);
  double maxInformation = shared_info_ptr_->maxInformationContent(goTermA);

  return sharedInformation / maxInformation;

}

