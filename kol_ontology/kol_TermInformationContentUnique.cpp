//
// Created by kellerberrin on 20/5/21.
//

#include "kol_TermInformationContentUnique.h"
#include "kel_exec_env.h"


namespace kol = kellerberrin::ontology;
namespace kel = kellerberrin;



void kol::TermInformationContentUnique::calcProbabilityMap( const std::shared_ptr<const GoGraph> &graph,
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
        annotations = annotations / const_bp_annotations;
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        annotations = annotations / const_mf_annotations;
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        annotations = annotations / const_cc_annotations;
        break;

      default:
      case GO::Ontology::ONTO_ERROR:
        ExecEnv::log().error("TermProbabilityUnique::calcProbabilityMap; GO term: {} does not have a valid ontology", term_id);
        break;

    }

  }

}
