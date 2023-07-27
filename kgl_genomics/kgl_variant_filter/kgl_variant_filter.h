///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant_db.h"
#include "kel_utility.h"
#include "kgl_variant_filter.h"
#include "kgl_variant_filter_type.h"


namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a specified minimum DP counts.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RefAltCountFilter : public FilterVariants {

public:

  explicit RefAltCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {

    std::stringstream ss;
    ss << "filter with minimum Ref+Alt base count:" << minimum_count_;
    filterName(ss.str());

  }
  ~RefAltCountFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<RefAltCountFilter>(*this); }

private:

  size_t minimum_count_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a specified minimum DP counts.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DPCountFilter : public FilterVariants {

public:

  explicit DPCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {

    std::stringstream ss;
    ss << "filter with minimum DP base count:" << minimum_count_;
    filterName(ss.str());

  }
  ~DPCountFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<DPCountFilter>(*this); }

private:

  size_t minimum_count_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants on phasing.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PhaseFilter : public FilterVariants {

public:

  explicit PhaseFilter(VariantPhase phase) : phase_(phase) {

    filterName("PhaseFilter");

  }
  ~PhaseFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.phaseId() == phase_; }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<PhaseFilter>(*this); }


private:

  const VariantPhase phase_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to SNPs (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PassFilter : public FilterVariants {

public:

  explicit PassFilter() { filterName("Variant marked 'Pass' for filters in VCF"); }
  ~PassFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.evidence().passFilter(); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<PassFilter>(*this); }


private:


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to SNPs (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPFilter : public FilterVariants {

public:

  explicit SNPFilter() { filterName("SNP and MNP filter"); }
  ~SNPFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.isSNP(); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<SNPFilter>(*this); }


private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// True Filter - performs no filtering. If combined with the NotFilter below, would filter every variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TrueFilter : public FilterVariants {

public:

  explicit TrueFilter() {

    filterName("TrueFilter");

  }
  ~TrueFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant&) const override { return true; }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<TrueFilter>(*this); }


private:

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// False Filter - unconditionally filters all variants. Useful for deleting large populations.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FalseFilter : public FilterVariants {

public:

  explicit FalseFilter() {

    filterName("FalseFilter");

  }
  ~FalseFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant&) const override { return false; }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FalseFilter>(*this); }

private:

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Negation Filter, the logical negation of a supplied filter.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class NotFilter : public FilterVariants {

public:

  explicit NotFilter(const FilterVariants& filter) {

    filter_ptr_ = std::dynamic_pointer_cast<FilterVariants>(filter.clone());
    filterName(std::string("NOT(") + filter_ptr_->filterName() + std::string(")"));

  }
  ~NotFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<NotFilter>(*this); }

private:


  std::shared_ptr<FilterVariants> filter_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// And Filter, logical and of two supplied filters.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AndFilter : public FilterVariants {

public:

  AndFilter(const FilterVariants& filter1, const FilterVariants& filter2) {

    filter1_ptr_ = std::dynamic_pointer_cast<FilterVariants>(filter1.clone());
    filter2_ptr_ = std::dynamic_pointer_cast<FilterVariants>(filter2.clone());

    filterName("AND(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")");

  }
  ~AndFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<AndFilter>(*this); }


private:

  std::shared_ptr<FilterVariants> filter1_ptr_;
  std::shared_ptr<FilterVariants> filter2_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Or Filter, logical or of two supplied filters.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class OrFilter : public FilterVariants {

public:

  OrFilter(const FilterVariants& filter1, const FilterVariants& filter2) {

    filter1_ptr_ = std::dynamic_pointer_cast<FilterVariants>(filter1.clone());
    filter2_ptr_ = std::dynamic_pointer_cast<FilterVariants>(filter2.clone());

    filterName("OR(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")");

  }
  ~OrFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<OrFilter>(*this); }

private:

  std::shared_ptr<FilterVariants> filter1_ptr_;
  std::shared_ptr<FilterVariants> filter2_ptr_;

};



}   // end namespace


#endif //KGL_FILTER_H
