//
// Created by kellerberrin on 24/04/23.
//

#ifndef KGL_VARIANT_FILTER_INFO_H
#define KGL_VARIANT_FILTER_INFO_H


#include "kgl_variant_db.h"
#include "kgl_variant_filter_db_variant.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"

#include <cctype>
#include <format>
#include <ranges>
#include <unordered_set>

namespace kellerberrin::genome {   //  organization::project level namespace


/// General Info filter class.
///
/// InfoType can only be templated with double, std::vector<double>, int64_t, std::vector<int64_t>, std::string,
/// std::vector<string> and bool.
/// Template Missing is the return value if the info field is not found.
template<typename InfoType, bool Missing>
requires ValidInfoDataType<InfoType>
class InfoFilter : public FilterVariants {

public:

  InfoFilter(std::string field_name, std::function<bool(const InfoType&)> filter_lambda)
      : field_name_(std::move(field_name)), filter_lambda_(std::move(filter_lambda)) {

    filterName(std::format("Info Filter: {}", field_name_));

  }
  ~InfoFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<InfoFilter>(*this); }

private:

  const std::string field_name_;
  const std::function<bool(const InfoType&)> filter_lambda_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename InfoType, bool Missing>
requires ValidInfoDataType<InfoType>
bool InfoFilter<InfoType, Missing>::applyFilter(const Variant& variant) const {

  auto info_opt = InfoEvidenceAnalysis::getTypedInfoData<InfoType>(variant, field_name_);
  if (info_opt) {

    return filter_lambda_(*info_opt);

  }

  return Missing;

}


/// Filter on a Vep subfield found in the Gnomad Homosapien data.
/// Note that if the variant has NO VEP evidence at all, then the filter returns false (NOT Missing).
/// If the vep field contains the specified sub-string then the filter returns 'true'.
/// If the empty string "" is specified then the corresponding vep field must be empty (all whitespace) to return true.
template<bool Missing>
class VepSubStringFilter : public FilterVariants {

public:

  VepSubStringFilter(std::string vep_field_name, std::string sub_string)
      : vep_field_name_(std::move(vep_field_name)),  sub_string_(std::move(sub_string)) {

    filterName(std::format("Vep Info SubField: {} contains sub string '{}'", vep_field_name_, sub_string_));

  }
  ~VepSubStringFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<VepSubStringFilter>(*this); }

private:

  const std::string vep_field_name_;
  const std::string sub_string_;

};


/// Vep Filter implementation.
template<bool Missing>
bool VepSubStringFilter<Missing>::applyFilter(const Variant& variant) const {

  auto vep_fields_opt = InfoEvidenceAnalysis::getVepSubFields(variant);

  if (not vep_fields_opt) {

    // No VEP evidence on this variant; intentionally treated as false, not Missing.
    return false;

  }

  const VEPSubFieldEvidence& vep_fields = *vep_fields_opt.value();

  auto vep_index_opt = vep_fields.vepHeader()->getSubFieldIndex(vep_field_name_);

  if (not vep_index_opt) {

    ExecEnv::log().error("VepSubStringFilter::applyFilter; could not find VEP field: {} in VEP fields", vep_field_name_);
    for (auto const& sub_field : vep_fields.vepHeader()->subFieldHeaders()) {

      ExecEnv::log().info("VepSubStringFilter::applyFilter; available VEP field: {} in VEP fields", sub_field);

    }
    return Missing;

  }

  size_t field_index = vep_index_opt.value();

  for (auto const& vep_field : vep_fields.vepFields()) {

    const std::vector<std::string_view>& sub_fields = VEPSubFieldEvidence::vepSubFields(vep_field);

    if (sub_fields.size() != vep_fields.vepHeader()->subFieldHeaders().size()) {

      ExecEnv::log().error("VepSubStringFilter::applyFilter; VEP sub-field count: {} not equal to VEP header size: {}",
                           sub_fields.size(), vep_fields.vepHeader()->subFieldHeaders().size());
      return Missing;

    }

    const std::string_view& sub_field = sub_fields[field_index];

    if (sub_string_.empty()) {

      // Equivalent to trimAllWhiteSpace(sub_field).empty() but without the per-field allocation.
      return std::ranges::all_of(sub_field, [](unsigned char c) { return std::isspace(c) != 0; });

    } else if (sub_field.find(sub_string_) != std::string::npos) {

      return true;

    }

  } // for all vep fields

  return false;

}


} // Namespace

#endif //KGL_VARIANT_FILTER_INFO_H
