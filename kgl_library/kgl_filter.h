///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant_vcf.h"
#include "kgl_genome_db.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a specified minimum DP counts.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RefAltCountFilter : public VariantFilter {

public:

  explicit RefAltCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {}
  ~RefAltCountFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const VCFVariant& variant) const override { return implementFilter(variant); }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<RefAltCountFilter>(*this); }

private:

  size_t minimum_count_;

  bool implementFilter(const Variant& variant) const;


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a specified minimum DP counts.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DPCountFilter : public VariantFilter {

public:

  explicit DPCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {}
  ~DPCountFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const VCFVariant& variant) const override { return implementFilter(variant); }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<DPCountFilter>(*this); }

private:

  size_t minimum_count_;

  bool implementFilter(const Variant& variant) const;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to SNPs (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPFilter : public VariantFilter {

public:

  explicit SNPFilter() {}
  ~SNPFilter() override = default;

  std::string filterName() const final { return "SNP and MNP Variants)"; }

  bool applyFilter(const VCFVariant& variant) const override { return variant.isSNP(); }

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

  bool applyFilter(const VCFVariant&) const override { return false; }

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

  bool applyFilter(const VCFVariant&) const override { return false; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InsertFilter>(*this); }

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFilter : public VariantFilter {

public:

  explicit ContigFilter(const ContigId_t& contig_ident) : contig_ident_(contig_ident) {}
  ~ContigFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const VCFVariant& variant) const override { return implementFilter(variant); }

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

  bool applyFilter(const VCFVariant& variant) const override { return variant.offset() >= start_ and variant.offset() < end_; }

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

  bool applyFilter(const VCFVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }

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

  bool applyFilter(const VCFVariant& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }

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

  bool applyFilter(const VCFVariant& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }

  std::string filterName() const final { return "OR(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")"; }

  std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<OrFilter>(*this); }

private:

  std::shared_ptr<VariantFilter> filter1_ptr_;
  std::shared_ptr<VariantFilter> filter2_ptr_;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_FILTER_H
