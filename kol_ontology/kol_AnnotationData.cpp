//
// Created by kellerberrin on 27/5/21.
//

#include "kel_exec_env.h"
#include "kol_AnnotationData.h"

namespace kol = kellerberrin::ontology;


void kol::AnnotationData::createAnnotationMap(const std::vector<std::shared_ptr<const GAFRecord>>& gaf_records) {

  for (auto const& gaf_record_ptr : gaf_records) {

    addGAFRecord(gaf_record_ptr);

  }

}


bool kol::AnnotationData::addGAFRecord(const std::shared_ptr<const GAFRecord>& gaf_record_ptr) {

  // Add term_id to map of gene_uniprot_id annotations.
  auto gene_iter = gene_annotation_map_.find(gaf_record_ptr->geneUniprotId());
  if (gene_iter == gene_annotation_map_.end()) {

    // Add the go record.
    GeneGOMap go_map;
    go_map.emplace(gaf_record_ptr->goIdent(), std::vector<std::shared_ptr<const GAFRecord>>{gaf_record_ptr});
    gene_annotation_map_.emplace(gaf_record_ptr->geneUniprotId(), std::move(go_map));

  } else {
    // Gene record exists.
    auto& [gene_id, go_map] = *gene_iter;
    // Check if go term exists for gene.
    auto go_iter = go_map.find(gaf_record_ptr->goIdent());
    if (go_iter == go_map.end()) {

      go_map.emplace(gaf_record_ptr->goIdent(), std::vector<std::shared_ptr<const GAFRecord>>{gaf_record_ptr});

    } else {

      auto& [go_term, gaf_record_vector] = *go_iter;
      gaf_record_vector.push_back(gaf_record_ptr);

    }

  }

  // Does go term exist.
  auto go_iter = go_annotation_map_.find(gaf_record_ptr->goIdent());
  if (go_iter == go_annotation_map_.end()) {

    // Add a gene record.
    GOGeneMap gene_map;
    gene_map.emplace(gaf_record_ptr->geneUniprotId(), std::vector<std::shared_ptr<const GAFRecord>>{gaf_record_ptr});
    go_annotation_map_.emplace(gaf_record_ptr->goIdent(), std::move(gene_map));

  } else {
    // Add a new gene record.
    auto& [go_id, gene_map] = *go_iter;
    auto gene_map_iter = gene_map.find(gaf_record_ptr->geneUniprotId());
    if (gene_map_iter == gene_map.end()) {

      gene_map.emplace(gaf_record_ptr->geneUniprotId(), std::vector<std::shared_ptr<const GAFRecord>>{gaf_record_ptr});

    } else {

      auto& [gene_id, gaf_record_vector] = *gene_map_iter;
      gaf_record_vector.push_back(gaf_record_ptr);

    }

  }

  return true;

}


std::vector<std::shared_ptr<const kol::GAFRecord>> kol::AnnotationData::filterGAFRecords(const PolicyEvidence& evidence_policy,
                                                                                            const std::vector<std::shared_ptr<const GAFRecord>>& gaf_records) {
  std::vector<std::shared_ptr<const GAFRecord>> filtered_records;

  for (auto const& gaf_record_ptr : gaf_records) {

    if (evidence_policy.isAllowed(gaf_record_ptr->evidenceCode())) {

      filtered_records.push_back(gaf_record_ptr);

    }

  }

  return filtered_records;

}


std::vector<std::string> kol::AnnotationData::getGoTermsForGene(const std::string &gene) const {

  std::vector<std::string> go_terms;
  for (auto const& [go_id, gaf_vector] : getGeneGOMap(gene)) {

    go_terms.push_back(go_id);

  }

  return go_terms;

}


std::vector<std::string> kol::AnnotationData::getGenesForGoTerm(const std::string &go_term) const {

  std::vector<std::string> genes;
  for (auto const& [gene_id, gaf_vector] : getGOGeneMap(go_term)) {

    genes.push_back(gene_id);

  }

  return genes;

}


std::vector<std::shared_ptr<const kol::GAFRecord>> kol::AnnotationData::getAllGAFRecords() const {

  std::vector<std::shared_ptr<const kol::GAFRecord>> gaf_records;

  for (auto const& [gene_id, go_map] : gene_annotation_map_) {

    for (auto const& [go_id, gaf_record_vector] : go_map) {

      for (auto const& gaf_record_ptr : gaf_record_vector) {

        gaf_records.push_back(gaf_record_ptr);

      }

    }

  }

  return gaf_records;

}




kol::OntologySetType<std::string> kol::AnnotationData::getGoTermsForGeneByOntology( const std::string &gene,
                                                                                       GO::Ontology filter_ontology,
                                                                                       const GoGraph &G) const {
  kol::OntologySetType<std::string> go_terms;

  for (auto const& [go_id, gaf_record_vector] : getGeneGOMap(gene)) {

    if (gaf_record_vector.empty()) {

      ExecEnv::log().error("AnnotationDataNew::getGoTermsForGeneByOntology; empty go record vector for gene id: {}, go id: {}", gene, go_id);
      return kol::OntologySetType<std::string>{};

    }
    if (gaf_record_vector.front()->ontology() == filter_ontology) {

      go_terms.insert(go_id);

    }

  }

  return go_terms;

}

const kol::GeneGOMap& kol::AnnotationData::getGeneGOMap(const std::string &gene) const {

  static const GeneGOMap empty_map;

  auto gene_iter = gene_annotation_map_.find(gene);
  if (gene_iter == gene_annotation_map_.end()) {

    return empty_map;

  } else {

    auto const& [gene_id, go_map] = *gene_iter;
    return go_map;

  }

}

const kol::GOGeneMap& kol::AnnotationData::getGOGeneMap(const std::string &go_ident) const {

  static const GOGeneMap empty_map;

  auto go_iter = go_annotation_map_.find(go_ident);
  if (go_iter == go_annotation_map_.end()) {

    return empty_map;

  } else {

    auto const& [go_id, gene_map] = *go_iter;
    return gene_map;

  }

}


std::vector<std::string> kol::AnnotationData::getOntologyTerms(const GoGraph &graph, GO::Ontology ontology) const {

  OntologySetType<std::string> go_term_set;

  for (auto const& [gene_id, go_map] : gene_annotation_map_) {

    for (auto const& [go_id, gaf_record_vector] : go_map) {

      if (gaf_record_vector.empty()) {

        ExecEnv::log().error("AnnotationDataNew::getOntologyTerms; empty go record vector for gene id: {}, go id: {}", gene_id, go_id);
        return std::vector<std::string>();

      }

      if (ontology == gaf_record_vector.front()->ontology()) {

        go_term_set.insert(go_id);

      }

    }

  }

  std::vector<std::string> go_term_vector;
  go_term_vector = SetUtilities::convertSet(go_term_set);

  return go_term_vector;

}

void kol::AnnotationData::addGenesForGoTerm(const std::string &go_ident, OntologySetType<std::string> &gene_set) const {


  for (auto const& [gene_id, gaf_vector] : getGOGeneMap(go_ident)) {

    gene_set.insert(gene_id);

  }

}
