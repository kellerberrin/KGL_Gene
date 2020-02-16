//
// Created by kellerberrin on 10/10/17.
//

#include "kel_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::Feature::addSubFeature(const FeatureIdent_t& sub_feature_id, std::shared_ptr<const Feature> sub_feature_ptr) {

  sub_features_.insert(std::make_pair(sub_feature_id, sub_feature_ptr));

}


void kgl::Feature::recusivelyPrintsubfeatures(long feature_level) const {

  for (const auto& feature : sub_features_) {

    ExecEnv::log().info("Level: {}, Feature: {}, Type: {}, begin: {}, end: {} strand: {}",
                        feature_level,
                        feature.second->id(),
                        feature.second->featureType(),
                        feature.second->sequence().begin(),
                        feature.second->sequence().end(),
                        static_cast<char>(feature.second->sequence().strand()));

    feature.second->recusivelyPrintsubfeatures(feature_level + 1); // Recursive call for the sub-feature.

  }

}



bool kgl::Feature::verifyCDSPhase(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) const {

  bool result = true;
  // Check for mod3
  for(const auto& sorted_cds : coding_seq_ptr->getMap()) {

    result = result and verifyMod3(sorted_cds.second->getSortedCDS());
    result = result and verifyStrand(sorted_cds.second->getSortedCDS());

  }

  return result;

}


bool kgl::Feature::verifyMod3(const SortedCDS& sorted_cds) const {

  bool result = true;
// Check the combined sequence length is mod 3 = 0

  ContigSize_t coding_sequence_length = 0;
  for (auto cds : sorted_cds) {

    coding_sequence_length += (cds.second->sequence().end() - cds.second->sequence().begin());

  }

  if ((coding_sequence_length % Codon::CODON_SIZE) != 0) {

    ExecEnv::log().warn("Gene: {} offset: {} CDS coding sequence length mod 3 not zero : {}",
                        id(),
                        sequence().begin(),
                        (coding_sequence_length % 3));

    result = false;

  }

  return result;

}

bool kgl::Feature::verifyStrand(const SortedCDS& sorted_cds) const {

  bool result = true;

// Check the strand is consistent and not unknown.
  for (auto cds : sorted_cds) {

    if (cds.second->sequence().strand() != sequence().strand()) {

      ExecEnv::log().error("CDS: {} offset: {} strand: {}, parent sequence strand: {} mis-match",
                           cds.second->id(),
                           cds.second->sequence().begin(),
                           static_cast<char>(cds.second->sequence().strand()),
                           static_cast<char>(sequence().strand()));
      result = false;

    }

  }

  return result;

}


std::shared_ptr<const kgl::Feature> kgl::Feature::getGene() const {

  // Recursively search upward

  if (hasSuperfeature()) {

    if (getSuperFeature()->isGene()) {

      return getSuperFeature();

    } else {

      return getSuperFeature()->getGene();

    }
  }

  return nullptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<const kgl::CodingSequenceArray>
kgl::GeneFeature::getCodingSequences(std::shared_ptr<const GeneFeature> gene) {

  std::shared_ptr<CodingSequenceArray> sequence_array_ptr(std::make_shared<CodingSequenceArray>());
  getCodingSequences(gene, gene, sequence_array_ptr);
  return sequence_array_ptr;

}

// This routine is recursive. Assumes all the CDS/EXONS are on the same sub-feature level.
bool kgl::GeneFeature::getCodingSequences(std::shared_ptr<const GeneFeature> gene_ptr,
                                          std::shared_ptr<const Feature> cds_parent_ptr,
                                          std::shared_ptr<CodingSequenceArray>& sequence_array_ptr) {

  bool result = true;
  SortedCDS parent_cds;

  for (auto sub_feature : cds_parent_ptr->subFeatures()) {

    if (sub_feature.second->isCDS()) {


      auto insert = parent_cds.insert(std::make_pair(sub_feature.second->sequence().begin(),
                                                     std::static_pointer_cast<const CDSFeature>(sub_feature.second)));

      if (not insert.second) {

        ExecEnv::log().warn("Duplicate coding feature: {} at contig offset: {}",
                            sub_feature.second->id(),
                            sub_feature.second->sequence().begin());
        result = false;
      }

    } else { // recursively call this function

      result = result and getCodingSequences(gene_ptr, sub_feature.second, sequence_array_ptr);

    }

  }

  if (not parent_cds.empty()) {

    std::shared_ptr<const CodingSequence> coding_sequence(std::make_shared<const CodingSequence>(gene_ptr,
                                                                                                 cds_parent_ptr,
                                                                                                 parent_cds));

    result = result and sequence_array_ptr->insertCodingSequence(coding_sequence);

  }

  return result;

}

