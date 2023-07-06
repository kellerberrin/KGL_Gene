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
// Filter candidate variants to only those that belong to any gene coding sequence defined by the Genome Reference.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FilterAllCodingVariants: public FilterVariants {

public:

  explicit FilterAllCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr) : all_coding_variants_(reference_ptr) {}
  ~FilterAllCodingVariants() override = default;
  FilterAllCodingVariants(const FilterAllCodingVariants& copy) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return all_coding_variants_.codingRegionVariant(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterAllCodingVariants>(*this); }

private:

  IntervalCodingVariants all_coding_variants_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter candidate variants to only those that belong to genes specified in the vector of GeneFeatures.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FilterCodingVariants: public FilterVariants {

public:

  explicit FilterCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector) : interval_coding_variants_(gene_vector) {}
  ~FilterCodingVariants() override = default;
  FilterCodingVariants(const FilterCodingVariants&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return interval_coding_variants_.codingRegionVariant(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterCodingVariants>(*this); }

private:

  IntervalCodingVariants interval_coding_variants_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter candidate variants to only those that belong to the specified GeneFeature.
// This includes 5 prime, coding, intron and 3 prime regions.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FilterGeneVariants: public FilterVariants {

public:

  explicit FilterGeneVariants(const std::shared_ptr<const GeneFeature>& gene_ptr) : gene_intervals_(gene_ptr) {}
  ~FilterGeneVariants() override = default;
  FilterGeneVariants(const FilterGeneVariants&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return gene_intervals_.isWithinGene(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterGeneVariants>(*this); }

private:

  GeneIntervalStructure gene_intervals_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter candidate variants to only those that belong to the CODING (CDS defined) transcripts of the specified GeneFeature.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FilterGeneCodingVariants: public FilterVariants {

public:

  explicit FilterGeneCodingVariants(const std::shared_ptr<const GeneFeature>& gene_ptr) : gene_intervals_(gene_ptr) {}
  ~FilterGeneCodingVariants() override = default;
  FilterGeneCodingVariants(const FilterGeneCodingVariants&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return gene_intervals_.codingModifier(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterGeneCodingVariants>(*this); }

private:

  GeneIntervalStructure gene_intervals_;

};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter candidate variants to only those that belong to a specific single CODING (CDS defined) transcript of the specified GeneFeature.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FilterGeneTranscriptVariants: public FilterVariants {

public:

  explicit FilterGeneTranscriptVariants(const std::shared_ptr<const GeneFeature>& gene_ptr, FeatureIdent_t transcript_id)
  : gene_intervals_(gene_ptr), transcript_id_(std::move(transcript_id)) {}
  ~FilterGeneTranscriptVariants() override = default;
  FilterGeneTranscriptVariants(const FilterGeneTranscriptVariants&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return gene_intervals_.transcriptModifier(variant, transcript_id_); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FilterGeneTranscriptVariants>(*this); }

private:

  GeneIntervalStructure gene_intervals_;
  FeatureIdent_t transcript_id_;

};


} // Namespace


#endif //KGL_VARIANT_FILTER_FEATURES_H
