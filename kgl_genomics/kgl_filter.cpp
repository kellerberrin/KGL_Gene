//
// Created by kellerberrin on 16/10/17.
//

#include <sstream>
#include <memory>
#include "kgl_filter.h"
#include "kel_utility.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"

namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Info Filter implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::InfoFilterImpl::applyFilter(const InfoFilterLambda& filter_lambda, const Variant& variant) const {

  auto field_data = InfoEvidenceAnalysis::getInfoData(variant, field_name_);

  if (field_data) {

    return filter_lambda(field_data.value());

  } else {

    return missing_default_;

  }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vep Filter implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::VepSubStringFilter::applyFilter(const Variant& variant) const {

  std::optional<std::unique_ptr<const VEPSubFieldEvidence>> vep_fields_opt = InfoEvidenceAnalysis::getVepSubFields(variant);

  if (vep_fields_opt) {

    const VEPSubFieldEvidence& vep_fields = *vep_fields_opt.value();

    std::optional<size_t> vep_index_opt = vep_fields.vepHeader()->getSubFieldIndex(vep_field_name_);

    if (not vep_index_opt) {

      ExecEnv::log().error("VepSubStringFilter::applyFilter; could not find VEP field: {} in VEP fields", vep_field_name_);
      for (auto const& sub_field : vep_fields.vepHeader()->subFieldHeaders()) {

        ExecEnv::log().info("VepSubStringFilter::applyFilter; available VEP field: {} in VEP fields", sub_field);

      }
      return false;

    }

    size_t field_index = vep_index_opt.value();

    for (auto const& vep_field : vep_fields.vepFields()) {

      const std::vector<std::string_view>& sub_fields = VEPSubFieldEvidence::vepSubFields(vep_field);

      if (sub_fields.size() != vep_fields.vepHeader()->subFieldHeaders().size()) {

        ExecEnv::log().error("VepSubStringFilter::applyFilter; VEP sub-field count: {} not equal to VEP header size: {}",
                             sub_fields.size(),vep_fields.vepHeader()->subFieldHeaders().size());
        return false;

      }

      const std::string_view& sub_field = sub_fields[field_index];

      if (sub_string_.empty()) {

        return kel::Utility::trimAllWhiteSpace(std::string(sub_field)).empty();

      } else if (sub_field.find(sub_string_) != std::string::npos) {

        return true;

      }

    } // for all vep fields

  } // has vep.

  return false;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a base count.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::RefAltCountFilter::implementFilter(const Variant& variant) const {

  std::optional<std::shared_ptr<const FormatData>> count_evidence_opt = variant.evidence().formatData();

  if (count_evidence_opt) {

    return (count_evidence_opt.value()->refCount() + count_evidence_opt.value()->altCount()) >= minimum_count_;

  } else {

    ExecEnv::log().info("RefAltCountFilter; variant does not have Ref+Alt base count evidence");

  }

  return true;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a base count.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::DPCountFilter::implementFilter(const Variant& variant) const {

  std::optional<std::shared_ptr<const FormatData>> count_evidence = variant.evidence().formatData();

  if (count_evidence) {

    return count_evidence.value()->DPCount() >= minimum_count_;

  } else {

    ExecEnv::log().info("DPCountFilter; variant does not have base count evidence");

  }

  return true;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::ContigFilter::implementFilter(const Variant& variant) const {

    return variant.contigId() == contig_ident_;

}




