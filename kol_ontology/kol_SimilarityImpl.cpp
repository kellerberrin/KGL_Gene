//
// Created by kellerberrin on 19/4/21.
//

#include "kol_SimilarityJiangConrath.h"
#include "kol_SimilarityLin.h"
#include "kol_SimilarityResnik.h"
#include "kol_SimilarityRelevance.h"

namespace kol = kellerberrin::ontology;


//! A method for calculating term-to-term similarity for GO terms using JiangConrath similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/

double kol::SimilarityJiangConrath::calculateTermSimilarityAlt(const std::string &go_termA, const std::string &go_termB) const {

  if (go_termA == go_termB) {

    return 1.0;

  }

  if (not shared_info_ptr_->validateTerms(go_termA, go_termB)) {

    return 0.0;

  }

  double sharedIC = shared_info_ptr_->sharedInformation(go_termA, go_termB);
  double termA_IC = shared_info_ptr_->termInformation(go_termA);
  double termB_IC = shared_info_ptr_->termInformation(go_termB);
  double maxIC = shared_info_ptr_->maxInformationContent(go_termA);

  double dist = termA_IC + termB_IC - (2.0 * sharedIC);

  return 1.0 - (dist / (2.0 * maxIC));

}

// This is the JiangConrath formula used by GOSemSim
double kol::SimilarityJiangConrath::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {

  //if the terms do not exist return 0.0 similarity
  if (not shared_info_ptr_->validateTerms(goTermA, goTermB)) {

    return 0.0;

  }

  double maxIC = shared_info_ptr_->maxInformationContent(goTermA);;
  double sharedIC = shared_info_ptr_->sharedInformation(goTermA, goTermB);;
  double termA_IC = shared_info_ptr_->termInformation(goTermA);
  double termB_IC = shared_info_ptr_->termInformation(goTermB);

  double dist = termA_IC + termB_IC - (2.0 * sharedIC);

  return 1.0 - std::min(1.0, dist/maxIC);;

}


//! A method for calculating term-to-term similarity for GO terms using Lin similarity
/*!
  This method returns the Resnik similarity or the information content of the most informative common ancestor.
*/
double kol::SimilarityLin::calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const {

  if (go_termA == go_termB) {

    return 1.0;

  }

  if (not shared_info_ptr_->validateTerms(go_termA, go_termB)) {

    return 0.0;

  }

  double sharedIC = shared_info_ptr_->sharedInformation(go_termA, go_termB);
  double termA_IC = shared_info_ptr_->termInformation(go_termA);
  double termB_IC = shared_info_ptr_->termInformation(go_termB);

  return (2.0 * sharedIC) / (termA_IC + termB_IC);

}


//! A method for calculating term-to-term similarity for GO terms using Normalized Resnik similarity
/*!
  This method returns the Resnik similarity divided by the maximum possible similarity
*/
double kol::SimilarityResnik::calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const {

  if (not shared_info_ptr_->validateTerms(go_termA, go_termB)) {

    return 0.0;

  }

  double shared_information = shared_info_ptr_->sharedInformation(go_termA, go_termB);
  double max_information = shared_info_ptr_->maxInformationContent(go_termA);

  return shared_information / max_information;

}


//! A method for calculating term-to-term similarity for GO terms using Relevance similarity
/*!
  This method returns the Relevance similarity.
*/
double kol::SimilarityRelevance::calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const {
  //if the terms do not exit return 0.0 similarity
  if (not shared_info_ptr_->validateTerms(go_termA, go_termB)) {

    return 0.0;

  }

  //get the MICA value (zero if no mica term)
  double mica_value = shared_info_ptr_->sharedInformation(go_termA, go_termB);
  double complement_prob_mica = 1.0 - std::exp(-1.0 * mica_value);
  double denom = (shared_info_ptr_->termInformation(go_termA) + shared_info_ptr_->termInformation(go_termB));

  if (denom == 0.0 or mica_value == 0.0) {

    return 0.0;

  }

  //return the normalized information content similarity of Relevance
  return ((2.0 * mica_value) / denom) * complement_prob_mica;

}

