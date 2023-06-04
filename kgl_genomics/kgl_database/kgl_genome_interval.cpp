//
// Created by kellerberrin on 29/05/23.
//

#include "kgl_genome_interval.h"



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>



namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::OpenRightInterval::OpenRightInterval(ContigOffset_t lower, ContigOffset_t upper)  {

  if (upper <= lower) {

    ExecEnv::log().warn("OpenRightInterval::OpenRightInterval, Incorrect Initialization, Upper Offset: {} <= Lower Offset: {}", upper, lower);
    if (upper == lower) {

      ++upper;

    } else {

      std::swap(lower, upper);

    }


  }

  lower_ = lower;
  upper_ = upper;

}


bool kgl::IntervalSet::containsInterval(const OpenRightInterval& interval) const {

  auto iter = this->lower_bound(interval);
  if (iter != this->end()) {

    if (iter->lower() == interval.lower()) {

      return iter->containsInterval(interval);

    }

  }

  iter = std::prev(iter, 1);
  if (iter != end()) {

    return iter->containsInterval(interval);

  }

  return false;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::GeneIntervalStructure::GeneIntervalStructure(const std::shared_ptr<const GeneFeature> &gene_feature) {

  gene_feature_ = gene_feature;
  gene_interval_ = OpenRightInterval(gene_feature->sequence().begin(), gene_feature->sequence().end());
  codingInterval(gene_feature);

}

void kgl::GeneIntervalStructure::codingInterval(const std::shared_ptr<const GeneFeature>& gene_feature) {

  // Get the gene transciptions.
  auto coding_sequence_array = GeneFeature::getTranscriptionSequences(gene_feature);

  for (auto const& [sequence_name, coding_sequence] : coding_sequence_array->getMap()) {

    IntervalSet sequence_intervals;
    for (auto const& [feature_offset, coding_feature] : coding_sequence->getFeatureMap()) {

      auto const& sequence = coding_feature->sequence();
      OpenRightInterval sequence_interval(sequence.begin(), sequence.end());

      // Add to the boost interval set.
      auto [insert_iter, result] = sequence_intervals.insert(sequence_interval);
      if (not result) {

        ExecEnv::log().warn("GeneIntervalStructure::codingInterval; Gene: {}, Transcript: {} has duplicate coding regions at offset: {}",
                            gene_feature->id(), sequence_name, sequence_interval.lower());

      }

    } // For each coding feature.

    auto [insert_iter, result] = gene_coding_transcripts_.try_emplace(sequence_name, sequence_intervals);
    if (not result) {

      ExecEnv::log().warn("GeneIntervalStructure::codingInterval; Gene: {}, unable to insert Transcript: {} (duplicate)",
                          gene_feature->id(), sequence_name);

    }

  } // For each coding sequence.

}


bool kgl::GeneIntervalStructure::isMemberCoding(ContigOffset_t offset) const {

  for (auto const& [transcript_id, interval_set] : gene_coding_transcripts_) {

    if (interval_set.containsOffset(offset)) {

      return true;

    }

  }

  return false;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::IntervalCodingVariants::IntervalCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector) {

  for (auto const& gene_feature : gene_vector) {

    auto contig_iter = contig_interval_map_.find(gene_feature->contig()->contigId());
    if (contig_iter == contig_interval_map_.end()) {

      auto [inserted_iter, result] = contig_interval_map_.try_emplace(gene_feature->contig()->contigId(), IntervalMap<GeneIntervalStructure>());
      if (not result) {

        ExecEnv::log().error("IntervalCodingVariants::IntervalCodingVariants; cannot insert contig: {} (duplicate)", gene_feature->contig()->contigId());
        continue;

      }

      contig_iter = inserted_iter;

    }

    auto& [contig_id, interval_map] = *contig_iter;

    GeneIntervalStructure gene_coding_intervals(gene_feature);
    auto gene_interval = gene_coding_intervals.geneInterval();
    auto [insert_iter, result] = interval_map.try_emplace(gene_interval, gene_coding_intervals);
    if (not result) {

      ExecEnv::log().error("IntervalCodingVariants::IntervalCodingVariants; cannot insert Gene : {}, Interval[{}, {}) (duplicate)"
                           , gene_feature->contig()->contigId(), gene_interval.lower(), gene_interval.upper());

    }

  } // For each gene.

}

bool kgl::IntervalCodingVariants::codingRegionVariant(const Variant &variant) const {

  // Check if the contig is defined.
  auto contig_iter = contig_interval_map_.find(variant.contigId());
  if (contig_iter == contig_interval_map_.end()) {

    return false;

  }

  auto const& [contig_id, interval_map] = *contig_iter;

  bool result{false};
  for (auto const& [gene_interval, interval_structure] : interval_map) {

    for (auto const& [transcript_id, coding_set] : interval_structure.codingTranscripts()) {

      // If the variant is not an delete indel then just check if it within a coding interval
      if (variant.variantType() != VariantType::INDEL_DELETE) {

        result = coding_set.containsOffset(variant.offset() + variant.alleleOffset());

      } else {
        // A delete indel and we need to check any upstream effect on a coding interval.
        ContigOffset_t delete_size = variant.reference().length() - variant.alternate().length();
        result = coding_set.containsOffset(variant.offset() + variant.alleleOffset() + delete_size);

      }

      if (result) {

        return true;

      }

    }

  }

  return result;

}


std::optional<std::shared_ptr<const kgl::GeneFeature>> kgl::IntervalCodingVariants::getGeneCoding(const Variant &variant) const {


  // Check if the contig is defined.
  auto contig_iter = contig_interval_map_.find(variant.contigId());
  if (contig_iter == contig_interval_map_.end()) {

    return std::nullopt;

  }

  auto const& [contig_id, interval_map] = *contig_iter;

  bool result{false};
  for (auto const& [gene_interval, interval_structure] : interval_map) {

    for (auto const& [transcript_id, coding_set] : interval_structure.codingTranscripts()) {

      // If the variant is not an delete indel then just check if it within a coding interval
      if (variant.variantType() != VariantType::INDEL_DELETE) {

        result = coding_set.containsOffset(variant.offset() + variant.alleleOffset());

      } else {
        // A delete indel and we need to check any upstream effect on a coding interval.
        ContigOffset_t delete_size = variant.reference().length() - variant.alternate().length();
        result = coding_set.containsOffset(variant.offset() + variant.alleleOffset() + delete_size);

      }

      if (result) {

        return interval_structure.getGene();

      }

    }

  }

  return std::nullopt;

}


std::optional<std::shared_ptr<const kgl::GeneFeature>> kgl::IntervalCodingVariants::getGeneInterval(const Variant &variant) const {


  // Check if the contig is defined.
  auto contig_iter = contig_interval_map_.find(variant.contigId());
  if (contig_iter == contig_interval_map_.end()) {

    return std::nullopt;

  }

  auto const& [contig_id, interval_map] = *contig_iter;

  // If the variant is not an delete indel then just check if it within a coding interval
  bool result{false};
  if (variant.variantType() != VariantType::INDEL_DELETE) {

    result = interval_map.containsOffset(variant.offset() + variant.alleleOffset());

  } else {
    // A delete indel and we need to check any upstream effect on a coding interval.
    ContigOffset_t delete_size = variant.reference().length() - variant.alternate().length();
    result = interval_map.containsOffset(variant.offset() + variant.alleleOffset() + delete_size);

  }

  if (result) {



  }

  return std::nullopt;

}

/*
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Pimpl implementation for filtering variants to arbitrary features.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::ImplementFeatureInterval::ImplementFeatureInterval(const std::vector<std::shared_ptr<const Feature>>& feature_vector) {


  for (auto const& feature : feature_vector) {

    auto contig_iter = contig_interval_map_.find(feature->contig()->contigId());
    if (contig_iter == contig_interval_map_.end()) {

      auto [inserted_iter, result] = contig_interval_map_.try_emplace(feature->contig()->contigId(), bil::interval_set<ContigOffset_t>());
      if (not result) {

        ExecEnv::log().error("FilterFeatureInterval::FilterFeatureInterval; cannot insert contig: {} (duplicate)", feature->contig()->contigId());
        continue;

      }

      contig_iter = inserted_iter;

    }

    auto& [contig_id, interval_set] = *contig_iter;

    auto const& sequence = feature->sequence();
    auto const sequence_interval{bil::discrete_interval<ContigOffset_t>::right_open(sequence.begin(), sequence.end())};

    // Add to the boost interval set.
    interval_set += sequence_interval;

  } // For each feature.

}

bool kgl::ImplementFeatureInterval::variantContained(const Variant &variant) const {

  // Check if the contig is defined.
  auto contig_iter = contig_interval_map_.find(variant.contigId());
  if (contig_iter == contig_interval_map_.end()) {

    return false;

  }

  auto const& [contig_id, interval_set] = *contig_iter;

  // For now just see if the variant is within the feature interval.
  bool result = bil::contains(interval_set, variant.offset() + variant.alleleOffset());

  return result;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Public facade classes.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::IntervalAllCodingVariants::IntervalAllCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr) {


  std::vector<std::shared_ptr<const GeneFeature>> gene_vector;
  for (auto const& [contig_id, contig_ptr] : reference_ptr->getMap()) {

    for (auto const& [gene_offset, gene_ptr] : contig_ptr->getGeneMap()) {

      gene_vector.push_back(gene_ptr);

    }

  }

  pimpl_coding_interval_ptr_ = std::make_shared<ImplementCodingInterval>(gene_vector);

}

kgl::IntervalAllCodingVariants::IntervalAllCodingVariants(const IntervalAllCodingVariants& copy) {

  pimpl_coding_interval_ptr_ = copy.pimpl_coding_interval_ptr_;

}

// Required due to forward decl of pimpl class.
kgl::IntervalAllCodingVariants::~IntervalAllCodingVariants() {}

bool kgl::IntervalAllCodingVariants::contains(const Variant& variant) const {

  return pimpl_coding_interval_ptr_->variantContained(variant);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::IntervalCodingVariants::IntervalCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector) {

  pimpl_coding_interval_ptr_ = std::make_shared<ImplementCodingInterval>(gene_vector);

}

kgl::IntervalCodingVariants::IntervalCodingVariants(const IntervalCodingVariants& copy) {

  pimpl_coding_interval_ptr_ = copy.pimpl_coding_interval_ptr_;

}

// Required due to forward decl of pimpl class.
kgl::IntervalCodingVariants::~IntervalCodingVariants() {}

bool kgl::IntervalCodingVariants::contains(const Variant& variant) const {

  return pimpl_coding_interval_ptr_->variantContained(variant);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::IntervalFeatures::IntervalFeatures(const std::vector<std::shared_ptr<const Feature>>& feature_vector) {

  pimpl_feature_interval_ptr_ = std::make_shared<ImplementFeatureInterval>(feature_vector);

}

kgl::IntervalFeatures::IntervalFeatures(const IntervalFeatures& copy) {

  pimpl_feature_interval_ptr_ = copy.pimpl_feature_interval_ptr_;

}

// Required due to forward decl of pimpl class.
kgl::IntervalFeatures::~IntervalFeatures() {}

bool kgl::IntervalFeatures::contains(const Variant& variant) const {

  return pimpl_feature_interval_ptr_->variantContained(variant);

}
*/
