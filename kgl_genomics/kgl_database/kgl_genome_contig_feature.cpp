//
// Created by kellerberrin on 1/12/18.
//


#include "kgl_genome_contig_feature.h"
#include "kgl_genome_contig.h"
#include "kel_patterns.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// StructuredFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::StructuredFeatures::addFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  id_feature_map_.insert(std::make_pair(feature_ptr->id(), feature_ptr));

  offset_feature_map_.insert(std::make_pair(feature_ptr->sequence().begin(), feature_ptr));

}

bool kgl::StructuredFeatures::findFeatureId(const FeatureIdent_t& feature_id,
                                           std::vector<std::shared_ptr<const Feature>>& feature_ptr_vec) const {

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

  for (auto [offset, feature_ptr] : offsetFeatureMap()) {

    Feature &feature = *feature_ptr;
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

      ExecEnv::log().warn("Feature: {}; has {} duplicate sub-features", feature.id(), duplicates);

    }

  }

}



void kgl::StructuredFeatures::removeSubFeatureDuplicates() {


  for (auto feature_pair : idFeatureMap()) {

    Feature &feature = *feature_pair.second;

    long duplicates_removed = deleteIterableDuplicates(feature.subFeatures());

    if (duplicates_removed > 0) {

      ExecEnv::log().info("{} duplicate sub-features removed from super feature: {}", duplicates_removed, feature.id());

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
  verifySubFeatureDuplicates();

}



bool kgl::StructuredFeatures::equivalent(const StructuredFeatures& lhs) const {

  if (offset_feature_map_.size() != lhs.offset_feature_map_.size()) {

    return false;

  }

  // Just compare the keys for now.
  auto lhs_offset_iter = lhs.offset_feature_map_.begin();
  for (auto const& [offset, offset_ptr] : offset_feature_map_) {

    if (lhs_offset_iter == lhs.offset_feature_map_.end()) {

      return false;

    }

    auto& [lhs_offset, lhs_offset_ptr] = *lhs_offset_iter;
    if (offset != lhs_offset and offset_ptr->equivalent(*lhs_offset_ptr)) {

      return false;

    }

    ++lhs_offset_iter;

  }

  if (id_feature_map_.size() != lhs.id_feature_map_.size()) {

    return false;

  }

  // Just compare the keys for now.
  auto lhs_id_iter = lhs.id_feature_map_.begin();
  for (auto const& [id, id_ptr] : id_feature_map_) {

    if (lhs_id_iter == lhs.id_feature_map_.end()) {

      return false;

    }

    auto& [lhs_id, lhs_id_ptr] = *lhs_id_iter;
    if (id != lhs_id and id_ptr->equivalent(*lhs_id_ptr)) {

      return false;

    }

    ++lhs_id_iter;

  }

  return true;

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

      std::vector<std::shared_ptr<const Feature>> super_feature_ptr_vec;

      if (not findFeatureId(super_feature_id, super_feature_ptr_vec)) {

        // Flag an Error; could not find super feature.
        ExecEnv::log().error("Feature: {}; Super Feature: {} does not exist", feature.id(), super_feature_id);

      }
      if (super_feature_ptr_vec.size() > 1) {

        // Warning, more than 1 super feature.
        ExecEnv::log().warn("Super Feature id: {} returned : {} Super Features",
                                 super_feature_id, super_feature_ptr_vec.size());

      }
      if (not super_feature_ptr_vec.empty()) {

        feature.setSuperFeature(super_feature_ptr_vec.front());
        std::const_pointer_cast<Feature>(super_feature_ptr_vec.front())->addSubFeature(feature.id(), feature_pair.second);

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

  std::vector<std::shared_ptr<const Feature>> feature_ptr_vec;

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

      const Feature &sub_feature = *sub_feature_pair.second;

      // TSS blocks can overlap superfeatures (Genes).
      if (sub_feature.isTSS()) continue;

      if (sub_feature.sequence().begin() < feature.sequence().begin()
          or sub_feature.sequence().end() > feature.sequence().end()) {

        ExecEnv::log().error("SubFeature: {}; [{}:{}] overlaps Feature {}; [{}:{}]",
                                  sub_feature.id(),
                                  sub_feature.sequence().begin(),
                                  sub_feature.sequence().end(),
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end());


      } // if sub-feature overlaps feature

    } // for all sub-features.

    // Check that features fit within the super-feature.
    if (feature.hasSuperfeature()) {

      const Feature &super_feature = *feature.getSuperFeature();

      if (feature.sequence().begin() < super_feature.sequence().begin()
          or feature.sequence().end() > super_feature.sequence().end()) {

        ExecEnv::log().error("Feature: {}; [{}:{}] overlaps SuperFeature {}; [{}:{}]",
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end(),
                                  super_feature.id(),
                                  super_feature.sequence().begin(),
                                  super_feature.sequence().end());

      } // If feature overlaps super-feature.

    }

  } // for all features.

}

