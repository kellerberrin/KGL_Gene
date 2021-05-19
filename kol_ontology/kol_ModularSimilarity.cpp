//
// Created by kellerberrin on 19/4/21.
//

#include "kol_ModularJiangConrath.h"
#include "kol_ModularLin.h"
#include "kol_ModularResnik.h"
#include "kol_ModularRelevance.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term-to-term similarity for GO terms using JiangConrath similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/

double kol::ModularJiangConrath::calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const {

  if (go_termA == go_termB) {

    return 1.0;

  }

  if (not shared_info_ptr_->validateTerms(go_termA, go_termB)) {

    return 0.0;

  }

  double sharedIC = shared_info_ptr_->sharedInformation(go_termA, go_termB);
  double termA_IC = shared_info_ptr_->sharedInformation(go_termA);
  double termB_IC = shared_info_ptr_->sharedInformation(go_termB);
  double maxIC = shared_info_ptr_->maxInformationContent(go_termA);

  double dist = termA_IC + termB_IC - (2.0 * sharedIC);

  return 1.0 - (dist / (2.0 * maxIC));

}


//! A method for calculating term-to-term similarity for GO terms using Lin similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/
double kol::ModularLin::calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const {

  if (go_termA == go_termB) {

    return 1.0;

  }

  if (not shared_info_ptr_->validateTerms(go_termA, go_termB)) {

    return 0.0;

  }

  double sharedIC = shared_info_ptr_->sharedInformation(go_termA, go_termB);
  double termA_IC = shared_info_ptr_->sharedInformation(go_termA);
  double termB_IC = shared_info_ptr_->sharedInformation(go_termB);

  return (2.0 * sharedIC) / (termA_IC + termB_IC);

}


//! A method for calculating term-to-term similarity for GO terms using Resnik similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/
double kol::ModularResnik::calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const {

  if (not shared_info_ptr_->validateTerms(go_termA, go_termB)) {

    return 0.0;

  }

  return shared_info_ptr_->sharedInformation(go_termA, go_termB);

}


//! A method for calculating term-to-term similarity for GO terms using Normalized Resnik similarity
/*!
  This method returns the Resnik similarity divided by the maximum possible similarity
*/
double kol::ModularResnik::calculateNormalizedTermSimilarity(const std::string &go_termA, const std::string &go_termB) const {

  if (not shared_info_ptr_->validateTerms(go_termA, go_termB)) {

    return 0.0;

  }

  double sharedInformation = shared_info_ptr_->sharedInformation(go_termA, go_termB);
  double maxInformation = shared_info_ptr_->maxInformationContent(go_termA);

  return sharedInformation / maxInformation;

}


//! A method for calculating term-to-term similarity for GO terms using Relevance similarity
/*!
  This method returns the Relevance similarity.
*/
double kol::ModularRelevance::calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const {
  //if the terms do not exit return 0.0 similarity
  if (not shared_info_ptr_->validateTerms(go_termA, go_termB)) {

    return 0.0;

  }

  //get the MICA value (zero if no mica term)
  double mica_value = shared_info_ptr_->sharedInformation(go_termA, go_termB);
  double complement_prob_mica = 1.0 - std::exp(-1.0 * mica_value);
  double denom = (shared_info_ptr_->sharedInformation(go_termA) + shared_info_ptr_->sharedInformation(go_termB));

  if (denom == 0.0 or mica_value == 0.0) {

    return 0.0;

  }

  //return the normalized information content similarity of Relevance
  return ((2.0 * mica_value) / denom) * complement_prob_mica;

}

