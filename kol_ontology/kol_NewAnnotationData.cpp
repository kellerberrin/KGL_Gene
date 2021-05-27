//
// Created by kellerberrin on 27/5/21.
//

#include "kel_exec_env.h"
#include "kol_NewAnnotationData.h"

namespace kol = kellerberrin::ontology;

bool kol::AnnotationDataNew::addAssociation( const std::string &gene_id,
                                             const std::string &go_term,
                                             GO::Ontology go_ontology,
                                             GO::EvidenceCode evidence_code) {



  // Add term_id to map of gene_id annotations.
  auto gene_iter = gene_annotation_map_.find(gene_id);
  if (gene_iter == gene_annotation_map_.end()) {

    // Add the go record.
    OntologyMapType<std::string, std::pair<GO::Ontology, GO::EvidenceCode>> go_map;
    go_map.emplace(go_term, std::pair<GO::Ontology, GO::EvidenceCode>{go_ontology, evidence_code});
    gene_annotation_map_.emplace(gene_id, std::move(go_map));

  } else {
    // Go record exists.
    auto& [gene_id, go_map] = *gene_iter;
    auto gene_go_iter = go_map.find(go_term);
    if (gene_go_iter == go_map.end()) {

      go_map.emplace(go_term, std::pair<GO::Ontology, GO::EvidenceCode>{go_ontology, evidence_code});

    } else {

      auto const& [go_term, ont_evidence_pair] = *gene_go_iter;
      auto const& [ontology, evidence] = ont_evidence_pair;
      if (go_ontology != ontology) {

        ExecEnv::log().error("AnnotationData::addAssociation; (Duplicate) Go term: {} to gene map: {}, with different ontology: {}"
                             , go_term, gene_id, GO::ontologyToString(ontology));
      }
      if (evidence_code != evidence) {

//        ExecEnv::log().error("AnnotationData::addAssociation; (Duplicate) Go term: {} to gene map: {}, with different evidence code: {}, existing GO evidence: {}"
//            , go_term, gene_id, GO::evidenceToString(evidence_code), GO::evidenceToString(evidence));
      }

    }

  }

  // Add gene_id to map of term_id annotations.
  auto go_iter = go_annotation_map_.find(go_term);
  if (go_iter == go_annotation_map_.end()) {

    // Add a gene record.
    OntologyMapType<std::string, GO::EvidenceCode> gene_map;
    gene_map.emplace(gene_id, evidence_code);
    std::pair<GO::Ontology, OntologyMapType<std::string, GO::EvidenceCode>> go_value{go_ontology, gene_map};
    go_annotation_map_.emplace(go_term, std::move(go_value));

  } else {
    // Add a new gene record.
    auto& [go_id, ont_map_pair] = *go_iter;
    auto& [ontology, gene_map] = ont_map_pair;
    auto [iter, result] = gene_map.try_emplace(gene_id, evidence_code);

    if (not result) {

//      ExecEnv::log().error("AnnotationData::addAssociation; problem adding (Duplicate) Gene: {} to go term map: {}", gene_id, go_term);

    }

  }

  return true;

}

