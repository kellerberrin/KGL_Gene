//
// Created by kellerberrin on 18/5/21.
//

#include "kol_OntologyTypes.h"
#include "kol_TermProbabilityUnique.h"
#include "kol_Accumulators.h"
#include "kel_exec_env.h"



namespace kol = kellerberrin::ontology;
namespace kel = kellerberrin;



double kol::TermProbabilityUnique::getRootCount(const std::string& root_id) const {

  auto const result = probability_map_.find(root_id);
  if (result == probability_map_.end()) {

    ExecEnv::log().error("TermProbabilityMap::getRootCount; root term: {} not in probability map", root_id);
    return 1.0;

  }

  auto const& [term_id, anno_ont_pair] = *result;
  auto const& [annotations, ontology] = anno_ont_pair;

  if (annotations <= 0.0) {

    ExecEnv::log().error("TermProbabilityMap::getRootCount; root term: {}, invalid annotation count: {}", root_id, annotations);
    return 1.0;

  }

  return annotations;

}



double kol::TermProbabilityUnique::getValue(const std::string &term_id) const {

  auto result = probability_map_.find(term_id);
  if (result == probability_map_.end()) {

    return badIdValue();

  }
  auto const& [map_term_id, prob_ont_pair] = *result;
  auto const& [probability, ontology] = prob_ont_pair;
  return probability;

}


//! Public method for calculating the most informative common ancestor value
/*!
  This method searches the sets to determine the most informative ancestor.
*/

double kol::TermProbabilityUnique::getMICAinfo(const OntologySetType<std::string> &ancestorsA,
                                               const OntologySetType<std::string> &ancestorsB) const {

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
double kol::TermProbabilityUnique::getEfficientMICA(const OntologySetType<std::string> &smaller_set,
                                                    const OntologySetType<std::string> &larger_set) const {

  double max{0.0};
  //loop over shorter list
  for (auto const &term : smaller_set) {

    if (larger_set.find(term) != larger_set.end()) {

      double term_value = getValue(term);
      if (term_value > max) {

        max = term_value;

      }

    }

  }

  return max;

}



void kol::TermProbabilityUnique::calcProbability(const std::shared_ptr<const GoGraph> &graph,
                                                 const std::shared_ptr<const AnnotationData> &annotation_data) {

  auto term_ont_map = graph->getAllOntTerms();
  for (auto const& [term_id, ontology] : term_ont_map) {

    size_t annotations{0};
    auto child_self_set = graph->getSelfDescendantTerms(term_id);
    for (auto const& child_term_id : child_self_set) {

      annotations += annotation_data->getNumAnnotationsForGoTerm(child_term_id);

    }

    std::pair<double, GO::Ontology> value_pair{static_cast<double>(annotations), ontology};
    auto [iter, result] = probability_map_.try_emplace(term_id, value_pair);
    if (not result) {

      ExecEnv::log().error("TermProbabilityUnique::calcProbability; Unable to add duplicate go term: {}", term_id);

    }

  }

  Accumulators::SimpleAccumulator minMaxBP;
  Accumulators::SimpleAccumulator minMaxMF;
  Accumulators::SimpleAccumulator minMaxCC;

  const double const_bp_annotations = getRootCount(GO::getRootTermBP());
  const double const_mf_annotations = getRootCount(GO::getRootTermMF());
  const double const_cc_annotations = getRootCount(GO::getRootTermCC());

  for (auto& [term_id, anno_ont_pair] : probability_map_) {

    auto& [annotations, ontology] = anno_ont_pair;

    if (annotations == 0.0) {

      continue;

    }

    switch (ontology) {

      case GO::Ontology::BIOLOGICAL_PROCESS:
        minMaxBP(annotations);
        annotations = annotations / const_bp_annotations;
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        minMaxMF(annotations);
        annotations = annotations / const_mf_annotations;
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        minMaxCC(annotations);
        annotations = annotations / const_cc_annotations;
        break;

      default:
      case GO::Ontology::ONTO_ERROR:
        ExecEnv::log().error("TermProbabilityUnique::calcProbabilityMap; GO term: {} does not have a valid ontology", term_id);
        break;

    }

  }

  //calculate single annotation minimum normalization factors
  bp_normalization_min_1anno_ = 1.0 / Accumulators::extractMax(minMaxBP);
  mf_normalization_min_1anno_ = 1.0 / Accumulators::extractMax(minMaxMF);
  cc_normalization_min_1anno_ = 1.0 / Accumulators::extractMax(minMaxCC);

  //calculate minimum annotation minimum normalization factors
  bp_normalization_min_min_anno_ = Accumulators::extractMin(minMaxBP) / Accumulators::extractMax(minMaxBP);
  mf_normalization_min_min_anno_ = Accumulators::extractMin(minMaxMF) / Accumulators::extractMax(minMaxMF);
  cc_normalization_min_min_anno_ = Accumulators::extractMin(minMaxCC) / Accumulators::extractMax(minMaxCC);

}
