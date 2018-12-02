//
// Created by kellerberrin on 1/12/18.
//


#include "kgl_genome_db.h"
#include "kgl_patterns.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// StructuredFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::StructuredFeatures::addFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  id_feature_map_.insert(std::make_pair(feature_ptr->id(), feature_ptr));

  offset_feature_map_.insert(std::make_pair(feature_ptr->sequence().begin(), feature_ptr));

}

bool kgl::StructuredFeatures::findFeatureId(const FeatureIdent_t& feature_id,
                                           std::vector<std::shared_ptr<kgl::Feature>>& feature_ptr_vec) const {

  auto iter_pair = id_feature_map_.equal_range(feature_id);

  feature_ptr_vec.clear();
  for (auto iter = iter_pair.first; iter != iter_pair.second; ++iter) {

    feature_ptr_vec.emplace_back(iter->second);

  }

  return not feature_ptr_vec.empty();

}


void kgl::StructuredFeatures::verifyContigOverlap() {

  // If feature dimensions are [1, size] instead of [0, size) then assume that conversion from the
  // Gff convention of [1, size] has not been performed correctly during feature read from disk.
  // Adjust to [0, size) here.
  // Note that this suggests a problem with the (3rd party) Gff read functionality and should be addressed there.

  for (auto feature_pair : offsetFeatureMap()) {

    Feature &feature = *feature_pair.second;
    // Error if feature overlaps the and of the contig.
    // If [1,contig_size] then adjust to [0, contig_size)

    if (feature.sequence().begin() == 1) { // adjust to [0, size)

      FeatureSequence adj_sequence = feature.sequence();
      adj_sequence.begin(0);
      feature.sequence(adj_sequence);
      ExecEnv::log().warn("Contig: {} 1-offset features [1, {}], adjusted to zero-offset [0, {})",
                          feature.contig()->contigId(), feature.contig()->contigSize(), feature.contig()->contigSize());

    } else if (feature.sequence().end() > feature.contig()->contigSize()) { // No features larger than the contig.

      FeatureSequence adj_sequence = feature.sequence();
      adj_sequence.end(feature.contig()->contigSize());
      feature.sequence(adj_sequence);
      ExecEnv::log().warn("Feature: {} [{}, {}) exceeds contig size :{} adjusted to [{}, {})",
                          feature.id(), feature.sequence().begin(), feature.sequence().end(),
                          feature.contig()->contigSize(), feature.sequence().begin(), feature.contig()->contigSize());

    }

  } // for contig

}


void kgl::StructuredFeatures::verifySubFeatureDuplicates() {

  for (auto feature_pair : idFeatureMap()) {
    Feature &feature = *feature_pair.second;

    long duplicates = checkDuplicates(feature.subFeatures());

    if (duplicates > 0) {

      kgl::ExecEnv::log().warn("Feature: {}; has {} duplicate sub-features", feature.id(), duplicates);

    }

  }

}


void kgl::StructuredFeatures::verifySuperFeatureDuplicates() {

  for (auto feature_pair : idFeatureMap()) {
    Feature &feature = *feature_pair.second;

    long duplicates = checkDuplicates(feature.superFeatures());

    if (duplicates > 0) {

      kgl::ExecEnv::log().warn("Feature: {}; has {} duplicate super-features", feature.id(), duplicates);

    }

  }

}



void kgl::StructuredFeatures::removeSubFeatureDuplicates() {


  for (auto feature_pair : idFeatureMap()) {

    Feature &feature = *feature_pair.second;

    long duplicates_removed = deleteIterableDuplicates(feature.subFeatures());

    if (duplicates_removed > 0) {

      kgl::ExecEnv::log().info("{} duplicate sub-features removed from super feature: {}", duplicates_removed, feature.id());

    }

  }

}

void kgl::StructuredFeatures::removeSuperFeatureDuplicates() {

  for (auto feature_pair : idFeatureMap()) {

    Feature &feature = *feature_pair.second;

    long duplicates_removed = deleteIterableDuplicates(feature.superFeatures());

    if (duplicates_removed > 0) {

      kgl::ExecEnv::log().info("{} duplicate super-features removed from sub-feature: {}", duplicates_removed, feature.id());

    }

  }

}


void kgl::StructuredFeatures::clearHierarchy() {

  // Remove all hierarchies for all features.
  for (auto feature_pair : idFeatureMap()) {

    Feature &feature = *feature_pair.second;
    // Remove feature hierarchy.

    feature.clearHierachy();

  }

}


void kgl::StructuredFeatures::verifyFeatureHierarchy() {

  verifyContigOverlap();
  removeSubFeatureDuplicates();
  removeSuperFeatureDuplicates();
  verifySubFeatureDuplicates();
  verifySuperFeatureDuplicates();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GeneExonFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::GeneExonFeatures::checkAddFeature(std::shared_ptr<Feature>& feature_ptr) {

  if (feature_ptr->isTSS()) {

    ExecEnv::log().error("GeneExonFeatures::checkAddFeature; Feature: {}, cannot add a TSS feature to the Gene feature structure",
                         feature_ptr->id());
    return false;

  }

  addFeature(feature_ptr);

  return true;

}


void kgl::GeneExonFeatures::setupVerifyHierarchy() {

  setupFeatureHierarchy();
  verifyGeneExonHierarchy();

}


void kgl::GeneExonFeatures::setupFeatureHierarchy() {

  // Remove all hierarchies for all features.
  clearHierarchy();

  // Establish or re-establish the hierarchies for all features.
  for (auto feature_pair : idFeatureMap()) {

    Feature& feature = *feature_pair.second;
    // For each feature lookup a list of super_features
    std::vector<FeatureIdent_t> super_features;

    feature.getAttributes().getSuperFeatureIds(super_features);

    // Add parent pointers for the child and child pointers for the super_features.
    for (auto super_feature_id : super_features) {

      std::vector<std::shared_ptr<kgl::Feature>> super_feature_ptr_vec;
      if (not findFeatureId(super_feature_id, super_feature_ptr_vec)) {

        // Flag an Error; could not find super feature.
        kgl::ExecEnv::log().error("Feature: {}; Super Feature: {} does not exist", feature.id(), super_feature_id);

      }
      if (super_feature_ptr_vec.size() > 1) {

        // Warning, more than 1 super feature.
        kgl::ExecEnv::log().warn("Super Feature id: {} returned : {} Super Features",
                                 super_feature_id, super_feature_ptr_vec.size());

      }
      for (auto& super_feature_ptr : super_feature_ptr_vec) {

        feature.addSuperFeature(super_feature_id, super_feature_ptr);
        super_feature_ptr->addSubFeature(feature.id(), feature_pair.second);

      } // For all super_features with same id.

    } // For all parent ids.

  } // For all features.

}



bool kgl::GeneExonFeatures::findGenes(ContigOffset_t offset, GeneVector &gene_ptr_vec) const {

  gene_ptr_vec.clear();
  auto lb_result = gene_map_.lower_bound(offset);

  if (lb_result == gene_map_.end()) {

    return false;

  }

  auto result = gene_map_.equal_range(lb_result->first);

  for (auto it = result.first; it != result.second; ++it) {

    if (offset >= it->second->sequence().begin()) {

      gene_ptr_vec.emplace_back(it->second);

    }

  }

  return not gene_ptr_vec.empty();

}



// Given a gene id and an mRNA id (sequence id) return the coding base sequence.
bool kgl::GeneExonFeatures::getCodingSequence(const FeatureIdent_t& gene_id,
                                            const FeatureIdent_t& sequence_id,
                                            std::shared_ptr<const CodingSequence>& coding_sequence_ptr) const {

  std::vector<std::shared_ptr<Feature>> feature_ptr_vec;
  std::shared_ptr<const GeneFeature> gene_ptr;
  if (findFeatureId(gene_id, feature_ptr_vec)) {

    for (const auto& feature_ptr : feature_ptr_vec) {

      if (feature_ptr->id() == gene_id) {

        if (feature_ptr->isGene()) {

          gene_ptr = std::dynamic_pointer_cast<const GeneFeature>(feature_ptr);
          break;

        } else {

          ExecEnv::log().warn("Feature: {} is not a gene.", gene_id);
          return false;

        }

      }

    }

  }

  if (not gene_ptr) {

    ExecEnv::log().warn("Gene not found for feature id: {}.", gene_id);
    return false;

  }

  std::shared_ptr<const CodingSequenceArray> sequence_array_ptr = GeneFeature::getCodingSequences(gene_ptr);

  if (sequence_array_ptr->empty()) {

    ExecEnv::log().warn("No valid coding sequences found for Gene: {}.", gene_id);
    return false;

  }

  for (const auto& sequence: sequence_array_ptr->getMap()) {

    if (sequence.second->getCDSParent()->id() == sequence_id) {

      coding_sequence_ptr = sequence.second;
      return true;

    }

  }

  ExecEnv::log().warn("No valid coding sequences found for sequence id: {}.", sequence_id);
  return false;

}


void kgl::GeneExonFeatures::verifyGeneExonHierarchy() {

  verifyFeatureHierarchy();
  verifySubFeatureSuperFeatureDimensions();
  createGeneMap();

}


void kgl::GeneExonFeatures::createGeneMap() {

  // Clear the lookup table.
  gene_map_.clear();

  // Iterate through all the features looking for Gene features.
  for(const auto& feature : offsetFeatureMap()) {

    if(feature.second->isGene()) {

      ContigOffset_t end_offset = feature.second->sequence().end();
      gene_map_.insert(std::make_pair(end_offset, std::static_pointer_cast<GeneFeature>(feature.second)));

    }

  }

}




void kgl::GeneExonFeatures::verifySubFeatureSuperFeatureDimensions() {


  // Check that sub-features fit within feature.
  for (auto feature_pair : idFeatureMap()) {

    Feature &feature = *feature_pair.second;

    // TSS blocks can overlap superfeatures (Genes).
    if (feature.isTSS()) continue;

    for (auto sub_feature_pair : feature.subFeatures()) {

      Feature &sub_feature = *sub_feature_pair.second;

      // TSS blocks can overlap superfeatures (Genes).
      if (sub_feature.isTSS()) continue;

      if (sub_feature.sequence().begin() < feature.sequence().begin()
          or sub_feature.sequence().end() > feature.sequence().end()) {

        kgl::ExecEnv::log().error("SubFeature: {}; [{}:{}] overlaps Feature {}; [{}:{}]",
                                  sub_feature.id(),
                                  sub_feature.sequence().begin(),
                                  sub_feature.sequence().end(),
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end());


      } // if sub-feature overlaps feature

    } // for all sub-features.

    // Check that features fit within super-feature.
    for (auto super_feature_pair : feature.superFeatures()) {
      Feature &super_feature = *super_feature_pair.second;
      if (feature.sequence().begin() < super_feature.sequence().begin()
          or feature.sequence().end() > super_feature.sequence().end()) {

        kgl::ExecEnv::log().error("Feature: {}; [{}:{}] overlaps SuperFeature {}; [{}:{}]",
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end(),
                                  super_feature.id(),
                                  super_feature.sequence().begin(),
                                  super_feature.sequence().end());

      } // If feature overlaps super-feature.

    } // for all super features.

  } // for all features.

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AdjalleyTSSFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::AdjalleyTSSFeatures::checkAddFeature(std::shared_ptr<Feature>& feature_ptr) {

  if (not feature_ptr->isTSS()) {

    ExecEnv::log().error("AdjalleyTSSFeatures::checkAddFeature; Feature: {}, can only add TSS features to the TSS feature structure",
                         feature_ptr->id());
    return false;

  }

  addFeature(feature_ptr);

  return true;

}


void kgl::AdjalleyTSSFeatures::setupVerifyHierarchy(const StructuredFeatures& gene_super_features) {

  setupFeatureHierarchy(gene_super_features);
  verifyFeatureHierarchy();

}


kgl::TSSVector kgl::AdjalleyTSSFeatures::getTSSVector() const {

  TSSVector tss_vector;

  for (auto feature : offsetFeatureMap()) {

    std::shared_ptr<const TSSFeature> tss_feature = std::dynamic_pointer_cast<TSSFeature>(feature.second);

    if (not tss_feature) {

      ExecEnv::log().error("Unexpected feature type for TSS feature: {}", feature.second->id());

    } else {

      tss_vector.push_back(tss_feature);

    }

  }

  return tss_vector;

}


void kgl::AdjalleyTSSFeatures::setupFeatureHierarchy(const StructuredFeatures& gene_super_features) {

  // Remove all hierarchies for all features.
  clearHierarchy();

  // Establish or re-establish the hierarchies for all features.
  for (auto feature_pair : idFeatureMap()) {

    Feature& feature = *feature_pair.second;

    // TSS features are assigned to GENES
    std::vector<FeatureIdent_t> assigned_features;
    feature.getAttributes().getAssignedFeatureIds(assigned_features);

    // Add parent pointers for the child and child pointers for the super_features.
    for (auto super_feature_id : assigned_features) {

      std::vector<std::shared_ptr<kgl::Feature>> super_feature_ptr_vec;
      if (not gene_super_features.findFeatureId(super_feature_id, super_feature_ptr_vec)) {

        // If the super feature is unassigned then continue.
        if (super_feature_id == TSSFeature::TSS_UNASSIGNED) {

          continue;

        }

        // Otherwise flag an Error; could not find super feature.
        kgl::ExecEnv::log().error("Feature: {}; Super Feature: {} does not exist", feature.id(), super_feature_id);

      }
      if (super_feature_ptr_vec.size() > 1) {

        // Warning, more than 1 super feature.
        kgl::ExecEnv::log().warn("Super Feature id: {} returned : {} Super Features",
                                 super_feature_id, super_feature_ptr_vec.size());

      }
      for (auto& super_feature_ptr : super_feature_ptr_vec) {

        feature.addSuperFeature(super_feature_id, super_feature_ptr);
        super_feature_ptr->addSubFeature(feature.id(), feature_pair.second);

      } // For all super_features with same id.

    } // For all parent ids.

  } // For all features.

}

