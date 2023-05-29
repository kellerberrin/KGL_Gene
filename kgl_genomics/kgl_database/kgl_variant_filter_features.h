//
// Created by kellerberrin on 15/05/23.
//

#ifndef KGL_VARIANT_FILTER_FEATURES_H
#define KGL_VARIANT_FILTER_FEATURES_H

#include "kgl_genome_genome.h"
#include "kgl_genome_interval.h"


namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter variants that are (possibly partially if an indel) within the coding sequence of a GENE.
// All genes within the reference genome are filtered.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FilterAllCodingVariants: public FilterVariants {

public:

  explicit FilterAllCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr) : all_coding_variants_(reference_ptr) {}
  ~FilterAllCodingVariants() override = default;

  FilterAllCodingVariants(const FilterAllCodingVariants&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return all_coding_variants_.contains(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterAllCodingVariants>(*this); }

private:

  IntervalAllCodingVariants all_coding_variants_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter variants that are (possibly partially if an indel) within the coding sequence of gene.
// The vector specifies which genes are filtered.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FilterCodingVariants: public FilterVariants {

public:

  explicit FilterCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector)
  : interval_coding_variants_(gene_vector) {}
  ~FilterCodingVariants() override = default;

  FilterCodingVariants(const FilterCodingVariants&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return interval_coding_variants_.contains(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterCodingVariants>(*this); }

private:

  IntervalCodingVariants interval_coding_variants_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter a vector of arbitrary features.
// If, for example, a feature is a gene, this can (and mostly does) include 5 prime, 3 prime and intron non-coding regions.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FilterFeatures: public FilterVariants {

public:

  explicit FilterFeatures(const std::vector<std::shared_ptr<const Feature>>& feature_vector)
  : interval_features_(feature_vector) {}
  ~FilterFeatures() override = default;

  FilterFeatures(const FilterFeatures&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return interval_features_.contains(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterFeatures>(*this); }

private:

  IntervalFeatures interval_features_;

};



} // Namespace


#endif //KGL_VARIANT_FILTER_FEATURES_H
