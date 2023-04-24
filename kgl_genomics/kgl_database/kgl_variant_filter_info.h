//
// Created by kellerberrin on 24/04/23.
//

#ifndef KGL_VARIANT_FILTER_INFO_H
#define KGL_VARIANT_FILTER_INFO_H



#include "kgl_variant.h"
#include "kel_utility.h"
#include "kgl_variant_filter.h"

#include <unordered_set>

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

  explicit InfoFilterImpl(std::string field_name) : field_name_(std::move(field_name)) {}
  InfoFilterImpl(const InfoFilterImpl&) = default;
  ~InfoFilterImpl() = default;

  [[nodiscard]] bool applyFilter(const InfoFilterLambda& filter_lambda, const Variant& variant) const;

  [[nodiscard]] std::string fieldName() const { return field_name_; }

private:

  const std::string field_name_;


};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// General Info filter class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoFilter : public FilterVariants {

public:

  InfoFilter(const std::string& field_name, InfoFilterLambda filter_lambda)
      : info_filter_(field_name), filter_lambda_(std::move(filter_lambda)) {

    filterName(field_name);

  }
  ~InfoFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<InfoFilter>(*this); }

private:

  const InfoFilterImpl info_filter_;
  const InfoFilterLambda filter_lambda_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter on a Vep subfield found in the Gnomad Homosapien data. Will silently return false for all other data.
//
// If the vep field contains the specified sub-string then the filter returns 'true'.
// If the the empty string "" is specified then the corresponding vep field must be empty to return true.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VepSubStringFilter : public FilterVariants {

public:

  VepSubStringFilter(std::string vep_field_name, std::string sub_string)  // How the filter responds if the data is missing.
      : vep_field_name_(std::move(vep_field_name)),  sub_string_(std::move(sub_string)) {

    std::stringstream ss;
    ss << "Vep Info SubField: " << vep_field_name_ << " contains sub string \"" << sub_string_ <<"\"";
    filterName(ss.str());

  }
  ~VepSubStringFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<VepSubStringFilter>(*this); }

private:

  const std::string vep_field_name_;
  const std::string sub_string_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Integer Info Filter. Returns true if a scalar integer and greater or equal to the comparison value.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoGEQIntegerFilter : public FilterVariants {

public:

  InfoGEQIntegerFilter(const std::string& field_name, int64_t comparison_value)  // How the filter responds if the data is missing.
      : info_filter_(field_name), comparison_value_(comparison_value) {

    std::stringstream ss;
    ss << "VCF Info Field: " << field_name << " >= " << comparison_value;
    filterName(ss.str());

  }
  ~InfoGEQIntegerFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<InfoGEQIntegerFilter>(*this); }


private:

  const InfoFilterImpl info_filter_;
  const int64_t comparison_value_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_integer_vector = std::get_if<std::vector<int64_t>>(&data_variant);
    if (p_integer_vector != nullptr) {

      if (p_integer_vector->size() == 1) {

        return (p_integer_vector->front() >= this->comparison_value_);

      } else {

        return false;

      }

    } else {

      return false;

    }

  };

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Float Info Filter. Returns true if a scalar float and greater or equal to the comparison value.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoGEQFloatFilter : public FilterVariants {

public:

  InfoGEQFloatFilter(const std::string& field_name, double comparison_value)
      : info_filter_(field_name), comparison_value_(comparison_value) {

    std::stringstream ss;
    ss << "VCF Info Field: " << field_name << " >= " << comparison_value;
    filterName(ss.str());

  }
  ~InfoGEQFloatFilter() override = default;


  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<InfoGEQFloatFilter>(*this); }

private:

  const InfoFilterImpl info_filter_;
  const double comparison_value_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_float_vector = std::get_if<std::vector<double>>(&data_variant);
    if (p_float_vector != nullptr) {

      if (p_float_vector->size() == 1) {

        return (p_float_vector->front() >= this->comparison_value_);

      } else {

        return false;

      }

    } else {

      return false;

    }

  };

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// String Info Filter. Returns true if a string Info field contains a specified substring (including itself).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoSubStringFilter : public FilterVariants {

public:

  InfoSubStringFilter(const std::string& field_name, std::string sub_string, bool case_sensitive = false)  // How the filter responds if the data is missing.
      : info_filter_(field_name), sub_string_(std::move(sub_string)), case_sensitive_(case_sensitive) {

    std::stringstream ss;
    ss << "VCF Info Field: " << field_name << " contains sub string \"" << sub_string <<"\"";
    filterName(ss.str());
    if (not case_sensitive_) {

      sub_string_ = Utility::toupper(sub_string_);

    }

  }
  ~InfoSubStringFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<InfoSubStringFilter>(*this); }

private:

  const InfoFilterImpl info_filter_;
  std::string sub_string_;
  const bool case_sensitive_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_string_vector = std::get_if<std::vector<std::string>>(&data_variant);
    if (p_string_vector != nullptr) {

      if (not p_string_vector->empty()) {

        for (auto const& value : *p_string_vector) {

          bool result;
          if (case_sensitive_) {

            result = value.find(sub_string_) != std::string::npos;

          } else {

            result = Utility::toupper(value).find(sub_string_) != std::string::npos;

          }

          if (result) {

            return true;

          }

        }

        return false;

      } else {

        return false;

      }

    } else {

      return false;

    }

  };

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Boolean Filter. Returns true if a boolean info field is true. Can be used with the Negation filter below.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoBooleanFilter : public FilterVariants {

public:

  explicit InfoBooleanFilter(const std::string& field_name)  // How the filter responds if the data is missing.
      : info_filter_(field_name) {

    std::stringstream ss;
    ss << "VCF Boolean Info Field: " << field_name;
    filterName(ss.str());

  }
  ~InfoBooleanFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<InfoBooleanFilter>(*this); }

private:

  const InfoFilterImpl info_filter_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_bool = std::get_if<bool>(&data_variant);
    if (p_bool != nullptr) {

      return *p_bool;

    } else {

      return false;

    }

  };

};



} // Namespace

#endif //KGL_VARIANT_FILTER_INFO_H
