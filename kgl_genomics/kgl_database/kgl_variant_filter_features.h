//
// Created by kellerberrin on 15/05/23.
//

#ifndef KGL_VARIANT_FILTER_FEATURES_H
#define KGL_VARIANT_FILTER_FEATURES_H

#include "kgl_variant_filter_type.h"
#include "kgl_genome_interval.h"


namespace kellerberrin::genome {   //  organization::project level namespace



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implemented using multimap.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


class MultiFilterAllCodingVariants: public FilterVariants {

public:

  explicit MultiFilterAllCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr) : all_coding_variants_(reference_ptr) {}
  ~MultiFilterAllCodingVariants() override = default;
  MultiFilterAllCodingVariants(const MultiFilterAllCodingVariants& copy) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return all_coding_variants_.codingRegionVariant(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<MultiFilterAllCodingVariants>(*this); }

private:

  IntervalCodingVariants all_coding_variants_;

};

class MultiFilterCodingVariants: public FilterVariants {

public:

  explicit MultiFilterCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector) : interval_coding_variants_(gene_vector) {}
  ~MultiFilterCodingVariants() override = default;
  MultiFilterCodingVariants(const MultiFilterCodingVariants&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return interval_coding_variants_.codingRegionVariant(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<MultiFilterCodingVariants>(*this); }

private:

  IntervalCodingVariants interval_coding_variants_;

};


} // Namespace


#endif //KGL_VARIANT_FILTER_FEATURES_H
