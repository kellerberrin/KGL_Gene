//
// Created by kellerberrin on 29/05/23.
//

#ifndef KGL_GENOME_INTERVAL_H
#define KGL_GENOME_INTERVAL_H


#include "kgl_variant_filter_type.h"
#include "kgl_genome_genome.h"


namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter variants that are (possibly partially if an indel) within the coding sequence of a GENE.
// All genes within the reference genome are filtered.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ImplementCodingInterval; // Pimpl implementation for determining variant membership of gene coding intervals.

class IntervalAllCodingVariants {

public:

  explicit IntervalAllCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr);
  ~IntervalAllCodingVariants();

  IntervalAllCodingVariants(const IntervalAllCodingVariants&);

  [[nodiscard]] bool contains(const Variant& variant) const;

private:

  std::shared_ptr<const ImplementCodingInterval> pimpl_coding_interval_ptr_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter variants that are (possibly partially if an indel) within the coding sequence of gene.
// The vector specifies which genes are filtered.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


class IntervalCodingVariants {

public:

  explicit IntervalCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector);
  ~IntervalCodingVariants();

  IntervalCodingVariants(const IntervalCodingVariants&);

  [[nodiscard]] bool contains(const Variant& variant) const;

private:

  std::shared_ptr<const ImplementCodingInterval> pimpl_coding_interval_ptr_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter a vector of arbitrary features.
// If, for example, a feature is a gene, this can (and mostly does) include 5 prime, 3 prime and intron non-coding regions.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ImplementFeatureInterval; // Pimpl implementation for arbitrary feature intervals only.

class IntervalFeatures {

public:

  explicit IntervalFeatures(const std::vector<std::shared_ptr<const Feature>>& feature_vector);
  ~IntervalFeatures();

  IntervalFeatures(const IntervalFeatures&);

  [[nodiscard]] bool contains(const Variant& variant) const;

private:

  std::shared_ptr<const ImplementFeatureInterval> pimpl_feature_interval_ptr_;

};



} // Namespace



#endif //KGL_GENOME_INTERVAL_H
