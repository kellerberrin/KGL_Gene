///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant_vcf.h"
#include "kgl_genome_db.h"
#include "kgl_variant_factory_vcf_evidence.h"

namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// General Info Filter.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using InfoFilterLambda = std::function<bool(const InfoDataVariant&)>;
class InfoFilter : public VariantFilter {

public:

  InfoFilter(const InfoSubscribedField& info_field, InfoFilterLambda filter_lambda, bool missing_default = false)  // How the filter responds if the data is missing.
   : info_field_(info_field), filter_lambda_(std::move(filter_lambda)), missing_default_(missing_default) {}
  InfoFilter(const InfoFilter&) = default;
  ~InfoFilter() override = default;

  [[nodiscard]] std::string filterName() const final { return "The Info Data Filter"; }

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return implementFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoFilter>(*this); }

private:

  InfoSubscribedField info_field_;
  InfoFilterLambda filter_lambda_;
  bool missing_default_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Integer Info Filter. Returns true if a scalar integer and greater or equal to the comparison value.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoGEQIntegerFilter : public VariantFilter {

public:

  InfoGEQIntegerFilter(const InfoSubscribedField& info_field, int64_t comparison_value, bool missing_default = false)  // How the filter responds if the data is missing.
  : info_field_(info_field), comparison_value_(comparison_value), missing_default_(missing_default) {}
  InfoGEQIntegerFilter(const InfoGEQIntegerFilter&) = default;
  ~InfoGEQIntegerFilter() override = default;

  [[nodiscard]] std::string filterName() const final { return "Info Data Greater or Equal Integer Filter"; }

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return implementFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoGEQIntegerFilter>(*this); }

private:

  InfoSubscribedField info_field_;
  int64_t comparison_value_;
  bool missing_default_;
  InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_integer_vector = std::get_if<std::vector<int64_t>>(&data_variant);
    if (p_integer_vector != nullptr) {

      if (p_integer_vector->size() == 1) {

        return (p_integer_vector->front() >= this->comparison_value_);

      } else {

        return missing_default_;

      }

    } else {

      return missing_default_;

    }

  };

  [[nodiscard]] bool implementFilter(const Variant& variant) const;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Float Info Filter. Returns true if a scalar float and greater or equal to the comparison value.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoGEQFloatFilter : public VariantFilter {

public:

  InfoGEQFloatFilter(const InfoSubscribedField& info_field, double comparison_value, bool missing_default = false)  // How the filter responds if the data is missing.
  : info_field_(info_field), comparison_value_(comparison_value), missing_default_(missing_default) {}
  InfoGEQFloatFilter(const InfoGEQFloatFilter&) = default;
  ~InfoGEQFloatFilter() override = default;

  [[nodiscard]] std::string filterName() const final { return "Info Data Greater or Equal Float Filter"; }

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return implementFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoGEQFloatFilter>(*this); }

private:

  InfoSubscribedField info_field_;
  double comparison_value_;
  bool missing_default_;
  InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_float_vector = std::get_if<std::vector<double>>(&data_variant);
    if (p_float_vector != nullptr) {

      if (p_float_vector->size() == 1) {

        return (p_float_vector->front() >= this->comparison_value_);

      } else {

        return missing_default_;

      }

    } else {

      return missing_default_;

    }

  };

  [[nodiscard]] bool implementFilter(const Variant& variant) const;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a specified minimum DP counts.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RefAltCountFilter : public VariantFilter {

public:

  explicit RefAltCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {}
  RefAltCountFilter(const RefAltCountFilter&) = default;
  ~RefAltCountFilter() override = default;

  [[nodiscard]] std::string filterName() const final;

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return implementFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<RefAltCountFilter>(*this); }

private:

  size_t minimum_count_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a specified minimum DP counts.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DPCountFilter : public VariantFilter {

public:

  explicit DPCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {}
  DPCountFilter(const DPCountFilter&) = default;
  ~DPCountFilter() override = default;

  [[nodiscard]] std::string filterName() const final;

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return implementFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<DPCountFilter>(*this); }

private:

  size_t minimum_count_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to SNPs (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPFilter : public VariantFilter {

public:

  explicit SNPFilter() = default;
  SNPFilter(const SNPFilter&) = default;
  ~SNPFilter() override = default;

  [[nodiscard]] std::string filterName() const final { return "SNP and MNP Variants)"; }

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return variant.isSNP(); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<SNPFilter>(*this); }

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFilter : public VariantFilter {

public:

  explicit ContigFilter(const ContigId_t& contig_ident) : contig_ident_(contig_ident) {}
  ContigFilter(const ContigFilter&) = default;
  ~ContigFilter() override = default;

  [[nodiscard]] std::string filterName() const final;

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return implementFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<ContigFilter>(*this); }

private:

  const ContigId_t contig_ident_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Region filter - generally used as AndFilter(ContigFilter("ContigName"), RegionFilter(Start, End))
// Uses the half open interval convention [Begin, End) with zero offset (the first contig element has zero offset).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RegionFilter : public VariantFilter {

public:

  explicit RegionFilter(ContigOffset_t start, ContigOffset_t end) : start_(start), end_(end) {}
  RegionFilter(const RegionFilter&) = default;
  ~RegionFilter() override = default;

  [[nodiscard]] std::string filterName() const final;

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return variant.offset() >= start_ and variant.offset() < end_; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<RegionFilter>(*this); }

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
  NotFilter(const NotFilter&) = default;
  ~NotFilter() override = default;

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }

  [[nodiscard]] std::string filterName() const final { return "NOT(" + filter_ptr_->filterName() + ")"; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<NotFilter>(*this); }

private:

  std::shared_ptr<VariantFilter> filter_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// And Filter
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AndFilter : public VariantFilter {

public:

  explicit AndFilter(const VariantFilter& filter1, const VariantFilter& filter2) : filter1_ptr_(filter1.clone()), filter2_ptr_(filter2.clone()) {}
  AndFilter(const AndFilter&) = default;
  ~AndFilter() override = default;

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }

  [[nodiscard]] std::string filterName() const final { return "AND(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")"; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<AndFilter>(*this); }

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
  OrFilter(const OrFilter&) = default;
  ~OrFilter() override = default;

  [[nodiscard]] bool applyFilter(const VCFVariant& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }

  [[nodiscard]] std::string filterName() const final { return "OR(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")"; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<OrFilter>(*this); }

private:

  std::shared_ptr<VariantFilter> filter1_ptr_;
  std::shared_ptr<VariantFilter> filter2_ptr_;

};



}   // end namespace


#endif //KGL_FILTER_H
