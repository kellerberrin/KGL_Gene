//
// Created by kellerberrin on 20/5/21.
//

#include "kol_InformationContentImpl.h"
#include "kel_exec_env.h"

namespace kol = kellerberrin::ontology;
namespace kel = kellerberrin;



double kol::InformationContentImpl::getRootCount(const std::string& root_id) const {

  auto const result = probability_map_.find(root_id);
  if (result == probability_map_.end()) {

    ExecEnv::log().error("TermProbabilityMap::getRootCount; root term: {} not in probfailure map", root_id);
    return 1.0;

  }

  auto const& [term_id, anno_ont_pair] = *result;
  auto const& [annotations, root_onto] = anno_ont_pair;

  if (annotations <= 0.0) {

    ExecEnv::log().error("TermProbabilityMap::getRootCount; root term: {}, invalid annotation count: {}", root_id, annotations);
    return 1.0;

  }

  return annotations;

}


double kol::InformationContentImpl::termInformation(const std::string &term_id) const {

  auto result = probability_map_.find(term_id);
  if (result == probability_map_.end()) {

    return BAD_INFO_VALUE_;

  }
  auto const& [map_term_id, prob_ont_pair] = *result;
  auto const& [probability, ontology] = prob_ont_pair;
  return probability;

}



/*!
  This function converts the probability terms of the probfailure map to information content.
*/
void kol::InformationContentImpl::convertProbtoIC() {

  for (auto& [term_id, prob_ont_pair] : probability_map_) {

    auto& [probability, ontology] = prob_ont_pair;

    if (probability <= 0.0) {

      probability = BAD_INFO_VALUE_;

    } else {

      probability = -1.0 * std::log(probability);
      switch(ontology) {

        case GO::Ontology::BIOLOGICAL_PROCESS:
          max_bp_information_ = std::max(probability, max_bp_information_);
          break;

        case GO::Ontology::MOLECULAR_FUNCTION:
          max_mf_information_ = std::max(probability, max_mf_information_);
          break;

        case GO::Ontology::CELLULAR_COMPONENT:
          max_cc_information_ = std::max(probability, max_cc_information_);
          break;

        default:
        case GO::Ontology::ONTO_ERROR:
          break;

      }

    }

  }


}

// Find if the information map contains two terms and they have the same ontology.
bool kol::InformationContentImpl::validateTerms(const std::string &id_termA, const std::string &id_termB) const {

  auto result_A = probability_map_.find(id_termA);
  if (result_A == probability_map_.end()) {

    return false;

  }

  auto result_B = probability_map_.find(id_termB);
  if (result_B == probability_map_.end()) {

    return false;

  }

  auto const& [termA, val_ont_A] = *result_A;
  auto const& [valA, ontA] = val_ont_A;

  auto const& [termB, val_ont_B] = *result_B;
  auto const& [valB, ontB] = val_ont_B;

  return ontA == ontB;

}


//! Public method for calculating the most informative common ancestor value
/*!
  This method searches the sets to determine the most informative ancestor.
*/

double kol::InformationContentImpl::sharedInformation(const std::string& go_termA, const std::string& go_termB, const GoGraphImpl &graph) const {

  //create 2 sets
  OntologySetType<std::string> ancestorsA = graph.getSelfAncestorTerms(go_termA);
  OntologySetType<std::string> ancestorsB = graph.getSelfAncestorTerms(go_termB);

  if (ancestorsA.empty() or ancestorsB.empty()) {

    return 0.0;

  }

  // Choose the smaller and larger set for maximum efficiency
  if (ancestorsA.size() < ancestorsB.size()) {

    return getEfficientMICA(ancestorsA, ancestorsB);

  } else {

    return getEfficientMICA(ancestorsB, ancestorsA);

  }

}

//! Private method for calculating the most informative common ancestor value
/*!
  This method searches the sets to determine the most informative ancestor.
*/
double kol::InformationContentImpl::getEfficientMICA(const OntologySetType<std::string> &smaller_set,
                                                     const OntologySetType<std::string> &larger_set) const {

  double max{0.0};
  //loop over shorter list
  for (auto const &term : smaller_set) {

    if (larger_set.find(term) != larger_set.end()) {

      double term_value = termInformation(term);
      if (term_value > max) {

        max = term_value;

      }

    }

  }

  return max;

}

double kol::InformationContentImpl::getMaxInformation(GO::Ontology ontology) const {

  switch(ontology) {

    case GO::Ontology::BIOLOGICAL_PROCESS:
      return max_bp_information_;

    case GO::Ontology::MOLECULAR_FUNCTION:
      return max_mf_information_;

    case GO::Ontology::CELLULAR_COMPONENT:
      return max_cc_information_;

    default:
    case GO::Ontology::ONTO_ERROR:
      return 0.0;

  }

}


double kol::InformationContentImpl::maxInformationContent(const std::string& term_id) const {

  auto result = probability_map_.find(term_id);
  if (result == probability_map_.end()) {

    return BAD_INFO_VALUE_;

  }

  auto const& [term, val_ont_pair] = *result;
  auto const& [value, ontology] = val_ont_pair;

  return getMaxInformation(ontology);

}

