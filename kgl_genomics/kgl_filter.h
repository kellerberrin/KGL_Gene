///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant.h"
#include "kgl_genome_db.h"

namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// These objects filter on Info field criteria specified in the supplied filter lambda.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using InfoFilterLambda = std::function<bool(const InfoDataVariant&)>;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper class finds the Info field (if it exists).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoFilterImpl  {

public:

  explicit InfoFilterImpl(std::string field_name, bool missing_default = false) : field_name_(std::move(field_name)),
                                                                                  missing_default_(missing_default) {}
  InfoFilterImpl(const InfoFilterImpl&) = default;
  ~InfoFilterImpl() = default;

  [[nodiscard]] bool applyFilter(const InfoFilterLambda& filter_lambda, const Variant& variant) const;

  [[nodiscard]] bool missingDefault() const { return missing_default_; }

  [[nodiscard]] std::string fieldName() const { return field_name_; }

private:

  const std::string field_name_;
  const bool missing_default_;


};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// General Info filter class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoFilter : public VariantFilter {

public:

  InfoFilter(const std::string& field_name, InfoFilterLambda filter_lambda, bool missing_default = false)
  : info_filter_(field_name, missing_default), filter_lambda_(std::move(filter_lambda)) {

    filterName(field_name);

  }
  InfoFilter(const InfoFilter&) = default;
  ~InfoFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoFilter>(*this); }


private:

  const InfoFilterImpl info_filter_;
  const InfoFilterLambda filter_lambda_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter on a Vep subfield found in the Gnomad Homo Sapien data. Will return return missing_default for all other data.
//
// Note that a variant may have several vep fields representing the different genomic structures the variant occupies.
// The filter is an 'or' against these fields.
// If any one (of several) of the specified vep fields meets the filter specification then the filter returns 'true'.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VepSubStringFilter : public VariantFilter {

public:

  VepSubStringFilter(std::string vep_field_name, std::string sub_string, bool missing_default = false)  // How the filter responds if the data is missing.
  : vep_field_name_(std::move(vep_field_name)),  sub_string_(std::move(sub_string)), missing_default_(missing_default) {

    std::stringstream ss;
    ss << "Vep Info SubField: " << vep_field_name_ << " contains sub string \"" << sub_string_ <<"\"";
    filterName(ss.str());

  }
  VepSubStringFilter(const VepSubStringFilter&) = default;
  ~VepSubStringFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<VepSubStringFilter>(*this); }

private:

  const std::string vep_field_name_;
  const std::string sub_string_;
  const bool missing_default_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Integer Info Filter. Returns true if a scalar integer and greater or equal to the comparison value.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoGEQIntegerFilter : public VariantFilter {

public:

  InfoGEQIntegerFilter(const std::string& field_name, int64_t comparison_value, bool missing_default = false)  // How the filter responds if the data is missing.
  : info_filter_(field_name, missing_default), comparison_value_(comparison_value) {

    std::stringstream ss;
    ss << "VCF Info Field: " << field_name << " >= " << comparison_value;
    filterName(ss.str());

  }
  InfoGEQIntegerFilter(const InfoGEQIntegerFilter&) = default;
  ~InfoGEQIntegerFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoGEQIntegerFilter>(*this); }

private:

  const InfoFilterImpl info_filter_;
  const int64_t comparison_value_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_integer_vector = std::get_if<std::vector<int64_t>>(&data_variant);
    if (p_integer_vector != nullptr) {

      if (p_integer_vector->size() == 1) {

        return (p_integer_vector->front() >= this->comparison_value_);

      } else {

        return info_filter_.missingDefault();

      }

    } else {

      return info_filter_.missingDefault();

    }

  };

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Float Info Filter. Returns true if a scalar float and greater or equal to the comparison value.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoGEQFloatFilter : public VariantFilter {

public:

  InfoGEQFloatFilter(const std::string& field_name, double comparison_value, bool missing_default = false)  // How the filter responds if the data is missing.
  : info_filter_(field_name, missing_default), comparison_value_(comparison_value) {

    std::stringstream ss;
    ss << "VCF Info Field: " << field_name << " >= " << comparison_value;
    filterName(ss.str());

  }
  InfoGEQFloatFilter(const InfoGEQFloatFilter&) = default;
  ~InfoGEQFloatFilter() override = default;


  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoGEQFloatFilter>(*this); }

private:

  const InfoFilterImpl info_filter_;
  const double comparison_value_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_float_vector = std::get_if<std::vector<double>>(&data_variant);
    if (p_float_vector != nullptr) {

      if (p_float_vector->size() == 1) {

        return (p_float_vector->front() >= this->comparison_value_);

      } else {

        return info_filter_.missingDefault();

      }

    } else {

      return info_filter_.missingDefault();

    }

  };

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// String Info Filter. Returns true if a string Info field contains a specified substring (including itself).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoSubStringFilter : public VariantFilter {

public:

  InfoSubStringFilter(const std::string& field_name, std::string sub_string, bool missing_default = false)  // How the filter responds if the data is missing.
  : info_filter_(field_name, missing_default), sub_string_(std::move(sub_string)) {

    std::stringstream ss;
    ss << "VCF Info Field: " << field_name << " contains sub string \"" << sub_string <<"\"";
    filterName(ss.str());

  }
  InfoSubStringFilter(const InfoSubStringFilter&) = default;
  ~InfoSubStringFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoSubStringFilter>(*this); }

private:

  const InfoFilterImpl info_filter_;
  const std::string sub_string_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_string_vector = std::get_if<std::vector<std::string>>(&data_variant);
    if (p_string_vector != nullptr) {

      if (p_string_vector->size() == 1) {

        return (p_string_vector->front().find(sub_string_) != std::string::npos);

      } else {

        return info_filter_.missingDefault();

      }

    } else {

      return info_filter_.missingDefault();

    }

  };

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Boolean Filter. Returns true if a boolean info field is true. Can be used with the Negation filter below.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoBooleanFilter : public VariantFilter {

public:

  explicit InfoBooleanFilter(const std::string& field_name, bool missing_default = false)  // How the filter responds if the data is missing.
  : info_filter_(field_name, missing_default) {

    std::stringstream ss;
    ss << "VCF Boolean Info Field: " << field_name;
    filterName(ss.str());

  }
  InfoBooleanFilter(const InfoBooleanFilter&) = default;
  ~InfoBooleanFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoBooleanFilter>(*this); }

private:

  const InfoFilterImpl info_filter_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_bool = std::get_if<bool>(&data_variant);
    if (p_bool != nullptr) {

      return *p_bool;

    } else {

      return info_filter_.missingDefault();

    }

  };

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a specified minimum DP counts.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RefAltCountFilter : public VariantFilter {

public:

  explicit RefAltCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {

    std::stringstream ss;
    ss << "Variants with minimum Ref+Alt base count:" << minimum_count_;
    filterName(ss.str());

  }
  RefAltCountFilter(const RefAltCountFilter&) = default;
  ~RefAltCountFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }

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

  explicit DPCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {

    std::stringstream ss;
    ss << "Variants with minimum DP base count:" << minimum_count_;
    filterName(ss.str());

  }
  DPCountFilter(const DPCountFilter&) = default;
  ~DPCountFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }

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

  explicit SNPFilter() {

    filterName("SNP and MNP Variants");

  }
  SNPFilter(const SNPFilter&) = default;
  ~SNPFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.isSNP(); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<SNPFilter>(*this); }

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFilter : public VariantFilter {

public:

  explicit ContigFilter(const ContigId_t& contig_ident) : contig_ident_(contig_ident) {

    filterName(std::string("Variants in Contig: ") + contig_ident_);

  }
  ContigFilter(const ContigFilter&) = default;
  ~ContigFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }

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

  explicit RegionFilter(ContigOffset_t start, ContigOffset_t end) : start_(start), end_(end) {

    std::stringstream ss;
    ss << "Variant in the half-interval [" << start_ << ", " << end_ << ")";
    filterName(ss.str());

  }
  RegionFilter(const RegionFilter&) = default;
  ~RegionFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.offset() >= start_ and variant.offset() < end_; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<RegionFilter>(*this); }

private:

  ContigOffset_t start_;
  ContigOffset_t end_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// True Filter - performs no filtering. If combined with the NotFilter below, would filter every variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TrueFilter : public VariantFilter {

public:

  explicit TrueFilter() {

    filterName("TrueFilter");

  }
  TrueFilter(const TrueFilter&) = default;
  ~TrueFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant&) const override { return true; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<TrueFilter>(*this); }

private:

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Negation Filter, the logical negation of a supplied filter.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class NotFilter : public VariantFilter {

public:

  explicit NotFilter(const VariantFilter& filter) : filter_ptr_(filter.clone()) {

    filterName(std::string("NOT(") + filter_ptr_->filterName() + std::string(")"));

  }
  NotFilter(const NotFilter&) = default;
  ~NotFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return not filter_ptr_->applyFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<NotFilter>(*this); }

private:

  std::shared_ptr<VariantFilter> filter_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// And Filter, logical and of two supplied filters.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AndFilter : public VariantFilter {

public:

  explicit AndFilter(const VariantFilter& filter1, const VariantFilter& filter2) : filter1_ptr_(filter1.clone()), filter2_ptr_(filter2.clone()) {

    filterName("AND(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")");

  }
  AndFilter(const AndFilter&) = default;
  ~AndFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<AndFilter>(*this); }

private:

  std::shared_ptr<VariantFilter> filter1_ptr_;
  std::shared_ptr<VariantFilter> filter2_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Or Filter, logical or of two supplied filters.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class OrFilter : public VariantFilter {

public:

  explicit OrFilter(const VariantFilter& filter1, const VariantFilter& filter2) : filter1_ptr_(filter1.clone()), filter2_ptr_(filter2.clone()) {

    filterName("OR(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")");

  }
  OrFilter(const OrFilter&) = default;
  ~OrFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<OrFilter>(*this); }

private:

  std::shared_ptr<VariantFilter> filter1_ptr_;
  std::shared_ptr<VariantFilter> filter2_ptr_;

};



}   // end namespace


#endif //KGL_FILTER_H
