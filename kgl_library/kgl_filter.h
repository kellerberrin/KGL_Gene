///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant_single.h"
#include "kgl_variant_compound.h"
#include "kgl_genome_db.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to SNPs (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPFilter : public VariantFilter {

public:

  explicit SNPFilter() {}
  ~SNPFilter() override = default;

  std::string filterName() const final { return "SNP and MNP Variants)"; }

  bool applyFilter(const SNPVariant&) const override { return true; }
  bool applyFilter(const DeleteVariant&) const override { return false; }
  bool applyFilter(const InsertVariant&) const override { return false; }
  bool applyFilter(const CompoundDelete&) const override { return false; }
  bool applyFilter(const CompoundInsert&) const override { return false; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<SNPFilter>(*this); }

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter to Delete variants (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DeleteFilter : public VariantFilter {

public:

  explicit DeleteFilter() {}
  ~DeleteFilter() override = default;

  std::string filterName() const final { return "Delete Variants"; }

  bool applyFilter(const SNPVariant&) const override { return false; }
  bool applyFilter(const DeleteVariant&) const override { return true; }
  bool applyFilter(const InsertVariant&) const override { return false; }
  bool applyFilter(const CompoundDelete&) const override { return true; }
  bool applyFilter(const CompoundInsert&) const override { return false; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<DeleteFilter>(*this); }

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter to Insert variants (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InsertFilter : public VariantFilter {

public:

  explicit InsertFilter() {}
  ~InsertFilter() override = default;

  std::string filterName() const final { return "Insert Variants"; }

  bool applyFilter(const SNPVariant&) const override { return false; }
  bool applyFilter(const DeleteVariant&) const override { return false; }
  bool applyFilter(const InsertVariant&) const override { return true; }
  bool applyFilter(const CompoundDelete&) const override { return false; }
  bool applyFilter(const CompoundInsert&) const override { return true; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InsertFilter>(*this); }

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants by quality (-10log10 {prob variant is in error})
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class QualityFilter : public VariantFilter {

public:

  explicit QualityFilter(Phred_t quality) : quality_(quality) {}
  ~QualityFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const SNPVariant& variant) const override { return variant.quality() >= quality_; }
  bool applyFilter(const DeleteVariant& variant) const override { return variant.quality() >= quality_; }
  bool applyFilter(const InsertVariant& variant) const override { return variant.quality() >= quality_; }
  bool applyFilter(const CompoundDelete& variant) const override { return variant.quality() >= quality_; }
  bool applyFilter(const CompoundInsert& variant) const override { return variant.quality() >= quality_; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<QualityFilter>(*this); }

private:

  const Phred_t quality_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFilter : public VariantFilter {

public:

  explicit ContigFilter(const ContigId_t& contig_ident) : contig_ident_(contig_ident) {}
  ~ContigFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return implementFilter(variant); }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<ContigFilter>(*this); }

private:

  const ContigId_t contig_ident_;

  bool implementFilter(const Variant& variant) const;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Region filter - generally used as AndFilter(ContigFilter("ContigName"), RegionFilter(Start, End))
// Uses the half open interval convention [Begin, End) with zero offset (the first contig element has zero offset).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RegionFilter : public VariantFilter {

public:

  explicit RegionFilter(ContigOffset_t start, ContigOffset_t end) : start_(start), end_(end) {}
  ~RegionFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const SNPVariant& variant) const override { return variant.offset() >= start_ and variant.offset() < end_; }
  bool applyFilter(const InsertVariant& variant) const override { return variant.offset() >= start_ and variant.offset() < end_; }
  bool applyFilter(const DeleteVariant& variant) const override { return variant.offset() >= start_ and variant.offset() < end_; }
  bool applyFilter(const CompoundDelete& variant) const override { return variant.offset() >= start_ and variant.offset() < end_; }
  bool applyFilter(const CompoundInsert& variant) const override { return variant.offset() >= start_ and variant.offset() < end_; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<RegionFilter>(*this); }

private:

  ContigOffset_t start_;
  ContigOffset_t end_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Negation Filter
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class NotFilter : public VariantFilter {

public:

  explicit NotFilter(const VariantFilter& filter) : filter_ptr_(filter.clone()) {}
  ~NotFilter() override = default;

  bool applyFilter(const SNPVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return not filter_ptr_->applyFilter(variant); }

  std::string filterName() const final { return "NOT(" + filter_ptr_->filterName() + ")"; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<NotFilter>(*this); }

private:

  std::shared_ptr<VariantFilter> filter_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// And Filter
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AndFilter : public VariantFilter {

public:

  explicit AndFilter(const VariantFilter& filter1, const VariantFilter& filter2) : filter1_ptr_(filter1.clone()), filter2_ptr_(filter2.clone()) {}
  ~AndFilter() override = default;

  bool applyFilter(const SNPVariant& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }

  std::string filterName() const final { return "AND(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")"; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<AndFilter>(*this); }

private:

  std::shared_ptr<VariantFilter> filter1_ptr_;
  std::shared_ptr<VariantFilter> filter2_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Or Filter
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class OrFilter : public VariantFilter {

public:

  explicit OrFilter(const VariantFilter& filter1, const VariantFilter& filter2) : filter1_ptr_(filter1.clone()), filter2_ptr_(filter2.clone()) {}
  ~OrFilter() override = default;

  bool applyFilter(const SNPVariant& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }

  std::string filterName() const final { return "OR(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")"; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<OrFilter>(*this); }

private:

  std::shared_ptr<VariantFilter> filter1_ptr_;
  std::shared_ptr<VariantFilter> filter2_ptr_;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_FILTER_H
