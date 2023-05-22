//
// Created by kellerberrin on 15/05/23.
//

#include "kgl_variant_filter_features.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <boost/icl/interval_set.hpp>
#include <boost/icl/discrete_interval.hpp>


namespace bt = boost;
namespace bil = boost::icl;
namespace kgl = kellerberrin::genome;

// Uses the boost interval container library.
using ContigIntervalMap = std::map<kgl::ContigId_t, bil::interval_set<kgl::ContigOffset_t>>;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Pimpl implementation for filtering variants to gene coding regions.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class kgl::GeneCodingInterval {

public:

  explicit GeneCodingInterval(const FilterGeneVector& gene_vector);
  ~GeneCodingInterval() = default;

  [[nodiscard]] bool variantContained(const Variant &variant) const;

private:

  ContigIntervalMap contig_interval_map_;

};


kgl::GeneCodingInterval::GeneCodingInterval(const FilterGeneVector& gene_vector) {

  for (auto const& gene_feature : gene_vector) {

    auto contig_iter = contig_interval_map_.find(gene_feature->contig()->contigId());
    if (contig_iter == contig_interval_map_.end()) {

      auto [inserted_iter, result] = contig_interval_map_.try_emplace(gene_feature->contig()->contigId(), bil::interval_set<ContigOffset_t>());
      if (not result) {

        ExecEnv::log().error("FilterFeatureInterval::FilterFeatureInterval; cannot insert contig: {} (duplicate)", gene_feature->contig()->contigId());
        continue;

      }

      contig_iter = inserted_iter;

    }

    auto& [contig_id, interval_set] = *contig_iter;

    // Get the gene transciptions.
    auto coding_sequence_array = GeneFeature::getTranscriptionSequences(gene_feature);

    for (auto const& [sequence_name, coding_sequence] : coding_sequence_array->getMap()) {

      for (auto const& [feature_offset, coding_feature] : coding_sequence->getFeatureMap()) {

        auto const& sequence = coding_feature->sequence();
        auto const sequence_interval{bil::discrete_interval<ContigOffset_t>::right_open(sequence.begin(), sequence.end())};

        // Add to the boost interval set.
        interval_set += sequence_interval;

      } // For each coding feature.

    } // For each coding sequence.

  } // For each gene.

}

bool kgl::GeneCodingInterval::variantContained(const Variant &variant) const {

  // Check if the contig is defined.
  auto contig_iter = contig_interval_map_.find(variant.contigId());
  if (contig_iter == contig_interval_map_.end()) {

    return false;

  }

  auto const& [contig_id, interval_set] = *contig_iter;

  bool result{false};
  // If the variant is not an delete indel then just check if it within a coding interval
  if (variant.variantType() != VariantType::INDEL_DELETE) {

    result = bil::contains(interval_set, variant.offset() + variant.alleleOffset());

  } else {
  // A delete indel and we need to check any upstream effect on a coding interval.
    ContigOffset_t delete_size = variant.reference().length() - variant.alternate().length();
    result = bil::contains(interval_set, variant.offset() + variant.alleleOffset() + delete_size);

  }

  return result;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Pimpl implementation for filtering variants to arbitrary features.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class kgl::FilterFeatureInterval {

public:

  explicit FilterFeatureInterval(const FilterFeatureVector& feature_vector);
  ~FilterFeatureInterval() = default;

  [[nodiscard]] bool variantContained(const Variant &variant) const;

private:

  ContigIntervalMap contig_interval_map_;

};


kgl::FilterFeatureInterval::FilterFeatureInterval(const FilterFeatureVector& feature_vector) {


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

bool kgl::FilterFeatureInterval::variantContained(const Variant &variant) const {

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
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::FilterAllCodingVariants::FilterAllCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr) {

  filterName("Filter Coding Variants");

  FilterGeneVector gene_vector;
  for (auto const& [contig_id, contig_ptr] : reference_ptr->getMap()) {

    for (auto const& [gene_offset, gene_ptr] : contig_ptr->getGeneMap()) {

        gene_vector.push_back(gene_ptr);

    }

  }

  pimpl_coding_interval_ptr_ = std::make_shared<GeneCodingInterval>(gene_vector);

}

kgl::FilterAllCodingVariants::FilterAllCodingVariants(const FilterAllCodingVariants& copy)  : FilterVariants(copy) {

  filterName(copy.filterName());
  pimpl_coding_interval_ptr_ = copy.pimpl_coding_interval_ptr_;

}

// Required due to forward decl of pimpl class.
kgl::FilterAllCodingVariants::~FilterAllCodingVariants() {}

bool kgl::FilterAllCodingVariants::applyFilter(const Variant& variant) const {

  return pimpl_coding_interval_ptr_->variantContained(variant);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::FilterCodingVariants::FilterCodingVariants(const FilterGeneVector& gene_vector) {

  filterName("Filter Coding Variants");
  pimpl_coding_interval_ptr_ = std::make_shared<GeneCodingInterval>(gene_vector);

}

kgl::FilterCodingVariants::FilterCodingVariants(const FilterCodingVariants& copy)  : FilterVariants(copy) {

  filterName(copy.filterName());
  pimpl_coding_interval_ptr_ = copy.pimpl_coding_interval_ptr_;

}

// Required due to forward decl of pimpl class.
kgl::FilterCodingVariants::~FilterCodingVariants() {}

bool kgl::FilterCodingVariants::applyFilter(const Variant& variant) const {

  return pimpl_coding_interval_ptr_->variantContained(variant);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::FilterFeatures::FilterFeatures(const FilterFeatureVector& feature_vector) {

  filterName("Filter Features");
  pimpl_feature_interval_ptr_ = std::make_shared<FilterFeatureInterval>(feature_vector);

}

kgl::FilterFeatures::FilterFeatures(const FilterFeatures& copy)  : FilterVariants(copy) {

  filterName(copy.filterName());
  pimpl_feature_interval_ptr_ = copy.pimpl_feature_interval_ptr_;

}

// Required due to forward decl of pimpl class.
kgl::FilterFeatures::~FilterFeatures() {}

bool kgl::FilterFeatures::applyFilter(const Variant& variant) const {

  return pimpl_feature_interval_ptr_->variantContained(variant);

}

