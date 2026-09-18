///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_VARIANT_FILTER_DB_VARIANT_H
#define KGL_VARIANT_FILTER_DB_VARIANT_H

#include "kgl_variant_db.h"
#include "kgl_variant_filter_type.h"

#include <format>


namespace kellerberrin::genome {   //  organization::project level namespace

/// Filter variants to a specified minimum Ref+Alt base count.
class RefAltCountFilter : public FilterVariants {

public:

  explicit RefAltCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {

    filterName(std::format("filter with minimum Ref+Alt base count:{}", minimum_count_));

  }
  ~RefAltCountFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<RefAltCountFilter>(*this); }

private:

  const size_t minimum_count_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;

};


/// Filter variants to a specified minimum DP base count.
class DPCountFilter : public FilterVariants {

public:

  explicit DPCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {

    filterName(std::format("filter with minimum DP base count:{}", minimum_count_));

  }
  ~DPCountFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<DPCountFilter>(*this); }

private:

  const size_t minimum_count_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;

};


/// Filter variants on phasing.
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


/// Filter variants to those marked 'Pass' for filters in the VCF.
class PassFilter : public FilterVariants {

public:

  explicit PassFilter() { filterName("Variant marked 'Pass' for filters in VCF"); }
  ~PassFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.evidence().passFilter(); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<PassFilter>(*this); }

};


/// Filter variants to SNPs (single and compound).
class SNPFilter : public FilterVariants {

public:

  explicit SNPFilter() { filterName("SNP filter"); }
  ~SNPFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.isSNP(); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<SNPFilter>(*this); }

};


/// Filter to indels that are not mod3 in size (frameshift).
class FrameShiftFilter : public FilterVariants {

public:

  explicit FrameShiftFilter() { filterName("Frame shift filter"); }
  ~FrameShiftFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FrameShiftFilter>(*this); }

private:

  [[nodiscard]] bool implementFilter(const Variant& variant) const;

};


/// True filter - performs no filtering. If combined with the NotFilter below, would filter every variant.
class TrueFilter : public FilterVariants {

public:

  explicit TrueFilter() {

    filterName("TrueFilter");

  }
  ~TrueFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant&) const override { return true; }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<TrueFilter>(*this); }

};


/// False filter - unconditionally filters all variants. Useful for deleting large populations.
class FalseFilter : public FilterVariants {

public:

  explicit FalseFilter() {

    filterName("FalseFilter");

  }
  ~FalseFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant&) const override { return false; }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<FalseFilter>(*this); }

};


/// Negation filter, the logical negation of a supplied filter.
class NotFilter : public FilterVariants {

public:

  explicit NotFilter(const FilterVariants& filter)
      : filter_ptr_(std::static_pointer_cast<FilterVariants>(filter.clone())) {

    filterName(std::format("NOT({})", filter_ptr_->filterName()));

  }
  ~NotFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<NotFilter>(*this); }

private:

  const std::shared_ptr<FilterVariants> filter_ptr_;

};


/// And filter, logical and of two supplied filters.
class AndFilter : public FilterVariants {

public:

  AndFilter(const FilterVariants& filter1, const FilterVariants& filter2)
      : filter1_ptr_(std::static_pointer_cast<FilterVariants>(filter1.clone())),
        filter2_ptr_(std::static_pointer_cast<FilterVariants>(filter2.clone())) {

    filterName(std::format("AND({}, {})", filter1_ptr_->filterName(), filter2_ptr_->filterName()));

  }
  ~AndFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<AndFilter>(*this); }

private:

  const std::shared_ptr<FilterVariants> filter1_ptr_;
  const std::shared_ptr<FilterVariants> filter2_ptr_;

};


/// Or filter, logical or of two supplied filters.
class OrFilter : public FilterVariants {

public:

  OrFilter(const FilterVariants& filter1, const FilterVariants& filter2)
      : filter1_ptr_(std::static_pointer_cast<FilterVariants>(filter1.clone())),
        filter2_ptr_(std::static_pointer_cast<FilterVariants>(filter2.clone())) {

    filterName(std::format("OR({}, {})", filter1_ptr_->filterName(), filter2_ptr_->filterName()));

  }
  ~OrFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<OrFilter>(*this); }

private:

  const std::shared_ptr<FilterVariants> filter1_ptr_;
  const std::shared_ptr<FilterVariants> filter2_ptr_;

};



}   // end namespace


#endif //KGL_VARIANT_FILTER_DB_VARIANT_H
