//
// Created by kellerberrin on 15/05/23.
//

#ifndef KGL_VARIANT_FILTER_FEATURES_H
#define KGL_VARIANT_FILTER_FEATURES_H

#include "kgl_variant_filter_type.h"
#include "kgl_genome_genome.h"



namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter variants that are (possibly partially if an indel) within the coding sequence of a GENE.
// All genes within the reference genome are filtered.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneCodingInterval; // Pimpl implementation for determining variant membership of gene coding intervals.

class FilterAllCodingVariants: public FilterVariants {

public:

  explicit FilterAllCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr);
  ~FilterAllCodingVariants() override;

  FilterAllCodingVariants(const FilterAllCodingVariants&);

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterAllCodingVariants>(*this); }

private:

  std::shared_ptr<const GeneCodingInterval> pimpl_coding_interval_ptr_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter variants that are (possibly partially if an indel) within the coding sequence of gene.
// The vector specifies which genes are filtered.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


using FilterGeneVector = std::vector<std::shared_ptr<const GeneFeature>>;
class FilterCodingVariants: public FilterVariants {

public:

  explicit FilterCodingVariants(const FilterGeneVector& gene_vector);
  ~FilterCodingVariants() override;

  FilterCodingVariants(const FilterCodingVariants&);

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterCodingVariants>(*this); }

private:

  std::shared_ptr<const GeneCodingInterval> pimpl_coding_interval_ptr_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter a vector of arbitrary features.
// If, for example, a feature is a gene, this can (and mostly does) include 5 prime, 3 prime and intron non-coding regions.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FilterFeatureInterval; // Pimpl implementation for arbitrary feature intervals only.

using FilterFeatureVector = std::vector<std::shared_ptr<const Feature>>;
class FilterFeatures: public FilterVariants {

public:

  explicit FilterFeatures(const FilterFeatureVector& feature_vector);
  ~FilterFeatures() override;

  FilterFeatures(const FilterFeatures&);

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterFeatures>(*this); }

private:

  std::shared_ptr<const FilterFeatureInterval> pimpl_feature_interval_ptr_;

};



} // Namespace


#endif //KGL_VARIANT_FILTER_FEATURES_H
