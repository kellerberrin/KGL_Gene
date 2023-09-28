//
// Created by kellerberrin on 1/12/18.
//


#include "kgl_genome_contig_feature.h"
#include "kgl_genome_contig.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// StructuredFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::StructuredFeatures::addFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  // Store by feature id and feature BEGIN offset.
  id_feature_map_.emplace(feature_ptr->id(), feature_ptr);
  offset_feature_map_.emplace(feature_ptr->sequence().begin(), feature_ptr);

}

std::vector<std::shared_ptr<const kgl::Feature>> kgl::StructuredFeatures::findFeatureId(const FeatureIdent_t& feature_id) const {

  std::vector<std::shared_ptr<const Feature>> feature_ptr_vec;

  auto const [lower_eq, upper_eq] = id_feature_map_.equal_range(feature_id);
  for (auto const& [feature_ident, feature_ptr] : std::ranges::subrange(lower_eq, upper_eq)) {

    feature_ptr_vec.emplace_back(feature_ptr);

  }

  return feature_ptr_vec;

}


void kgl::StructuredFeatures::verifyContigOverlap() const {

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


size_t kgl::StructuredFeatures::verifySubFeatureDuplicates() const {

  size_t total_duplicates{0};

  for (auto const& [feature_ident, feature_ptr] : idFeatureMap()) {

    std::set<FeatureIdent_t> check_duplicate_set;
    for (auto const& [sub_feature_ident, sub_feature_ptr] : feature_ptr->subFeatures()) {

         check_duplicate_set.insert(sub_feature_ident);

    }

    size_t duplicate_count = feature_ptr->subFeatures().size() - check_duplicate_set.size();
    total_duplicates += duplicate_count;

    if (duplicate_count > 0) {

      ExecEnv::log().warn("StructuredFeatures::verifySubFeatureDuplicates; Feature: {}; has {} duplicate sub-features",
                          feature_ident, duplicate_count);

    }

  }

  return total_duplicates;

}


void kgl::StructuredFeatures::removeSubFeatureDuplicates() {

  size_t feature_duplicates = verifySubFeatureDuplicates();

  if (feature_duplicates > 0) {

    for (auto const& [feature_ident, feature_ptr] : idFeatureMap()) {

      std::map<std::string, std::shared_ptr<const Feature>> unique_features_map;
      for (auto const& [sub_feature_ident, sub_feature_ptr] : feature_ptr->subFeatures()) {

        unique_features_map[sub_feature_ident] = sub_feature_ptr;

      }

      size_t duplicates_removed = feature_ptr->subFeatures().size() - unique_features_map.size();
      if (duplicates_removed > 0) {

        ExecEnv::log().info("StructuredFeatures::removeSubFeatureDuplicates; {} duplicate sub-features removed from super feature: {}",
                            duplicates_removed, feature_ident);

        feature_ptr->subFeatures().clear();
        for (auto const& [sub_feature_ident, sub_feature_ptr] : unique_features_map) {

          feature_ptr->subFeatures().emplace(sub_feature_ident, sub_feature_ptr);

        }

      } // If sub-feature duplicates

    }  // For all sub-features.

  } // If any duplicates

}


void kgl::StructuredFeatures::clearHierarchy() {

  // Remove all hierarchies for all features.
  for (auto const& [feature_ident, feature_ptr] : idFeatureMap()) {

    // Remove feature hierarchy.
    feature_ptr->clearHierachy();

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
  for (auto const& [feature_ident, feature_ptr] : idFeatureMap()) {

    // For each feature lookup a list of super_features
    std::vector<FeatureIdent_t> super_features;
    feature_ptr->getAttributes().getSuperFeatureIds(super_features);

    // Add parent pointers for the child and child pointers for the super_features.
    for (auto const& super_feature_id : super_features) {

      std::vector<std::shared_ptr<const Feature>> super_ptr_vec = findFeatureId(super_feature_id);
      if (super_ptr_vec.empty()) {

        // Flag an Error; could not find super feature.
        ExecEnv::log().error("GeneExonFeatures::setupFeatureHierarchy; Feature: {}; Super Feature: {} does not exist",
                             feature_ident, super_feature_id);

      }
      if (super_ptr_vec.size() > 1) {

        // Warning, more than 1 super feature.
        ExecEnv::log().warn("GeneExonFeatures::setupFeatureHierarchy; Super Feature id: {} returned : {} Super Features",
                            super_feature_id, super_ptr_vec.size());

      }
      if (not super_ptr_vec.empty()) {

        feature_ptr->setSuperFeature(super_ptr_vec.front());
        std::const_pointer_cast<Feature>(super_ptr_vec.front())->addSubFeature(feature_ident, feature_ptr);

      } // For all super_features with same id.

    } // For all parent ids.

  } // For all features.

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
  for(const auto& [offset, feature_ptr] : offsetFeatureMap()) {

    if (feature_ptr->isGene()) {

      // Note that genes are indexed by their END offset.
      ContigOffset_t end_offset = feature_ptr->sequence().end();
      gene_map_.emplace(end_offset, std::static_pointer_cast<GeneFeature>(feature_ptr));

    }

  }

}



void kgl::GeneExonFeatures::verifySubFeatureSuperFeatureDimensions() {


  // Check that sub-features fit within feature.
  for (auto const& [feature_id, feature_ptr] : idFeatureMap()) {

    for (auto const& [sub_feature_id, sub_feature_ptr] : feature_ptr->subFeatures()) {

      if (not checkSubFeatures(feature_ptr, sub_feature_ptr)) {

        ExecEnv::log().warn("SubFeature: {}, type: {}; {}[{}:{}) overlaps Feature {} type: {}; {}[{}:{})",
                                  sub_feature_ptr->id(),
                            sub_feature_ptr->type(),
                                  sub_feature_ptr->sequence().strandText(),
                                  sub_feature_ptr->sequence().begin(),
                                  sub_feature_ptr->sequence().end(),
                                  feature_ptr->id(),
                            feature_ptr->type(),
                                  feature_ptr->sequence().strandText(),
                                  feature_ptr->sequence().begin(),
                                  feature_ptr->sequence().end());


      } // if sub-feature overlaps feature

    } // for all sub-features.

    // Check that features fit within the super-feature.
    if (feature_ptr->hasSuperfeature()) {

      auto super_feature_ptr = feature_ptr->getSuperFeature();

      if (not checkSuperFeature(feature_ptr, super_feature_ptr)) {

        ExecEnv::log().warn("Feature: {}, type: {}; {}[{}:{}) overlaps SuperFeature {}, type: {}; {}[{}:{})",
                                  feature_ptr->id(),
                            feature_ptr->type(),
                                  feature_ptr->sequence().strandText(),
                                  feature_ptr->sequence().begin(),
                                  feature_ptr->sequence().end(),
                                  super_feature_ptr->id(),
                            super_feature_ptr->type(),
                                  super_feature_ptr->sequence().strandText(),
                                  super_feature_ptr->sequence().begin(),
                                  super_feature_ptr->sequence().end());

      } // If feature overlaps super-feature.

    }

  } // for all features.

}


bool kgl::GeneExonFeatures::checkSubFeatures( const std::shared_ptr<const Feature>& feature_ptr
                                             , const std::shared_ptr<const Feature>& sub_feature_ptr) {

  // TSS blocks can overlap features.
  if (sub_feature_ptr->isTSS()) {

    return true;

  }

  // 5' UTR should be at the beginning of the mRNA feature but can overlap.
  if (sub_feature_ptr->isUTR5() and feature_ptr->ismRNA()) {

    return true;

  }

  // 3' UTR should be at the end of the mRNA feature but can overlap.
  if (sub_feature_ptr->isUTR3() and feature_ptr->ismRNA()) {

    return true;

  }

  // Else just check the sub_feature is within the feature.
  bool valid_subfeature = sub_feature_ptr->sequence().begin() >= feature_ptr->sequence().begin()
  and sub_feature_ptr->sequence().end() <= feature_ptr->sequence().end();

  return valid_subfeature;

}


bool kgl::GeneExonFeatures::checkSuperFeature( const std::shared_ptr<const Feature>& feature_ptr
    , const std::shared_ptr<const Feature>& super_feature_ptr) {

  // TSS blocks can overlap superfeatures.
  if (feature_ptr->isTSS()) {

    return true;

  }

  // 5' UTR should be at the beginning of the mRNA feature but can overlap.
  if (feature_ptr->isUTR5() and super_feature_ptr->ismRNA()) {

    return true;

  }

  // 3' UTR should be at the end of the mRNA feature but can overlap.
  if (feature_ptr->isUTR3() and super_feature_ptr->ismRNA()) {

    return true;

  }

  // Else just check the feature is within the superfeature.
  bool valid_feature = feature_ptr->sequence().begin() >= super_feature_ptr->sequence().begin()
      and feature_ptr->sequence().end() <= super_feature_ptr->sequence().end();

  return valid_feature;

}