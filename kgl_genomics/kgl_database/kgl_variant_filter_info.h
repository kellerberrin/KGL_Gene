//
// Created by kellerberrin on 24/04/23.
//

#ifndef KGL_VARIANT_FILTER_INFO_H
#define KGL_VARIANT_FILTER_INFO_H



#include "kgl_variant.h"
#include "kel_utility.h"
#include "kgl_variant_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"

#include <unordered_set>

namespace kellerberrin::genome {   //  organization::project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// General Info filter class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// InfoType can only be templated with double, std::vector<double>, int64_t, std::vector<int64_t>, std::string,
// std::vector<string> and bool.
// Template Missing is the return value if the info field is not found.

template<typename InfoType, bool Missing>
requires ValidInfoDataType<InfoType>
class InfoFilter : public FilterVariants {

public:

  InfoFilter(const std::string& field_name, const std::function<bool(const InfoType&)>& filter_lambda)
      : field_name_(field_name), filter_lambda_(filter_lambda) {

    filterName("Info Filter: " + field_name);

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

  auto info_opt = InfoEvidenceAnalysis::getTypedInfoData<InfoType>( variant, field_name_);
  if (info_opt) {

    return filter_lambda_(info_opt.value());

  } else {

    return Missing;

  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter on a Vep subfield found in the Gnomad Homosapien data. Will silently return false for all other data.
//
// If the vep field contains the specified sub-string then the filter returns 'true'.
// If the empty string "" is specified then the corresponding vep field must be empty to return true.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<bool Missing>
class VepSubStringFilter : public FilterVariants {

public:

  VepSubStringFilter(std::string vep_field_name, std::string sub_string)
      : vep_field_name_(std::move(vep_field_name)),  sub_string_(std::move(sub_string)) {

    std::stringstream ss;
    ss << "Vep Info SubField: " << vep_field_name_ << " contains sub string '" << sub_string_ << "'";
    filterName(ss.str());

  }
  ~VepSubStringFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<VepSubStringFilter>(*this); }

private:

  const std::string vep_field_name_;
  const std::string sub_string_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vep Filter implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<bool Missing>
bool VepSubStringFilter<Missing>::applyFilter(const Variant& variant) const {

  std::optional<std::unique_ptr<const VEPSubFieldEvidence>> vep_fields_opt = InfoEvidenceAnalysis::getVepSubFields(variant);

  if (vep_fields_opt) {

    const VEPSubFieldEvidence& vep_fields = *vep_fields_opt.value();

    std::optional<size_t> vep_index_opt = vep_fields.vepHeader()->getSubFieldIndex(vep_field_name_);

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
                             sub_fields.size(),vep_fields.vepHeader()->subFieldHeaders().size());
        return Missing;

      }

      const std::string_view& sub_field = sub_fields[field_index];

      if (sub_string_.empty()) {

        return Utility::trimAllWhiteSpace(std::string(sub_field)).empty();

      } else if (sub_field.find(sub_string_) != std::string::npos) {

        return true;

      }

    } // for all vep fields

  } // has vep.

  return false;

}


} // Namespace

#endif //KGL_VARIANT_FILTER_INFO_H
