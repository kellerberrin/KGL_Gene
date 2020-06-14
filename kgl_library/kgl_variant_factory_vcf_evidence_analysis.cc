//
// Created by kellerberrin on 5/6/20.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"

#include <vector>

namespace kgl = kellerberrin::genome;


std::vector<std::string> kgl::InfoEvidenceAnalysis::varianttoStrings(const InfoDataVariant& info_data) {

  auto p_string_vector = std::get_if<std::vector<std::string>>(&info_data);
  if (p_string_vector != nullptr) {

    return *p_string_vector;

  } else {

    ExecEnv::log().error("InfoEvidenceAnalysis::varianttoStrings, Info data variant does not contain a string vector.");
    return std::vector<std::string>();

  }

}


std::vector<double> kgl::InfoEvidenceAnalysis::varianttoFloats(const InfoDataVariant& info_data) {

  auto p_float_vector = std::get_if<std::vector<double>>(&info_data);
  if (p_float_vector != nullptr) {

    return *p_float_vector;

  } else {

    ExecEnv::log().error("InfoEvidenceAnalysis::varianttoFloats, Info data variant does not contain a float vector.");
    return std::vector<double>();

  }

}



std::vector<int64_t> kgl::InfoEvidenceAnalysis::varianttoIntegers(const InfoDataVariant& info_data) {

  auto p_integer_vector = std::get_if<std::vector<int64_t>>(&info_data);
  if (p_integer_vector != nullptr) {

    return *p_integer_vector;

  } else {

    ExecEnv::log().error("InfoEvidenceAnalysis::varianttoIntegers, Info data variant does not contain a integer vector.");
    return std::vector<int64_t>();

  }

}


bool kgl::InfoEvidenceAnalysis::variantToBool(const InfoDataVariant& info_data) {

  auto p_bool = std::get_if<bool>(&info_data);
  if (p_bool != nullptr) {

    return *p_bool;

  } else {

    ExecEnv::log().error("InfoEvidenceAnalysis::varianttoBool, Info data variant does not contain a boolean.");
    return false;

  }

}


// Converts a bin in string format "1|0|0|0|1|0|0|0|1|0" into a vector of floats.
// If an error then returns a zero vector of expected size.
std::vector<double> kgl::InfoEvidenceAnalysis::stringBinToFloat(const std::vector<std::string>& bin_data, size_t expected_bin_size) {

  if (bin_data.size() >= 1) {

    std::vector<std::string> bin_strings = Utility::char_tokenizer(bin_data.front(), BIN_DELIMITER_);

    if (bin_strings.size() != expected_bin_size) {

      ExecEnv::log().warn("InfoEvidenceAnalysis::stringBinToFloat, Expected Bin Size: {}, Actual Bin Size: {}, Bin String",
                           expected_bin_size, bin_strings.size(), bin_data.front());
      return std::vector<double>(expected_bin_size, 0.0);

    }

    std::vector<double> bin_vector;

    for (auto const& bin : bin_strings) {

      try {

        bin_vector.emplace_back(std::stod(bin));

      }
      catch (...) {

        ExecEnv::log().error( "InfoEvidenceAnalysis::stringBinToFloat, problem converting bin: {} to double, bin string: {}",
                              bin, bin_data.front());
        return std::vector<double>(expected_bin_size, 0.0);

      }

    }

    return bin_vector;

  } else {

    return std::vector<double>(expected_bin_size, 0.0);

  }

}


std::optional<const kgl::InfoSubscribedField> kgl::InfoEvidenceAnalysis::getSubscribedField( const Variant& variant,
                                                                                             const std::string& field_ident) {

  if (variant.evidence().infoData()) {

    return variant.evidence().infoData().value()->evidenceHeader()->getSubscribedField(field_ident);

  }

  return std::nullopt;

}


std::optional<kgl::InfoDataVariant> kgl::InfoEvidenceAnalysis::getInfoData( const Variant& variant,
                                                                            const std::string& field_ident) {

  std::optional<const kgl::InfoSubscribedField> field = getSubscribedField(variant, field_ident);

  if (field) {

    const DataMemoryBlock &data_block = *variant.evidence().infoData().value();

    InfoDataVariant variant_data = field->getData(data_block);

    return variant_data;

  }

  return std::nullopt;

}

// Creates a vep subfield object, if defined by the variant.
std::optional<std::unique_ptr<const kgl::VEPSubFieldEvidence>>
kgl::InfoEvidenceAnalysis::getVepSubFields(const Variant& variant) {

  std::optional<const kgl::InfoSubscribedField> vep_field_opt = getSubscribedField( variant, VEPSubFieldHeader::VEP_FIELD_ID);

  if (not vep_field_opt) {

    return std::nullopt;

  }

  std::optional<std::shared_ptr<const VEPSubFieldHeader>> vep_header_opt = vep_field_opt.value().vepSubFieldHeader();

  if (not vep_header_opt) {

    ExecEnv::log().error("InfoEvidenceAnalysis::getVepSubFields, found 'vep' field but no matching 'vep' header defined");
    return std::nullopt;

  }

  std::shared_ptr<const VEPSubFieldHeader> vep_header_ptr = vep_header_opt.value();

  if (vep_header_ptr->subFieldHeaders().empty()) {

    ExecEnv::log().error("InfoEvidenceAnalysis::getVepSubFields, found empty 'vep' header");
    return std::nullopt;

  }

  std::optional<kgl::InfoDataVariant> vep_opt = getInfoData(variant, VEPSubFieldHeader::VEP_FIELD_ID);

  if (not vep_opt) {

    return std::nullopt;

  }

  std::vector<std::string> vep_field_vector = varianttoStrings(vep_opt.value());

  // Only add the vep fields that correspond to the variant alternate.
  // There may be multiple alternates in the vep vector.
  std::vector<std::string> variant_vep_vector;
  variant_vep_vector.reserve(vep_field_vector.size());
  std::vector<std::vector<std::string_view>> vep_sub_fields_vector;
  vep_sub_fields_vector.reserve(vep_field_vector.size());
  for (auto& vep_field : vep_field_vector) {

    std::vector<std::string_view> vep_sub_fields = Utility::view_tokenizer(vep_field, VEPSubFieldHeader::VEP_DELIMITER_CHAR);
    if (vep_header_ptr->subFieldHeaders().size() != vep_sub_fields.size()) {

      ExecEnv::log().error( "InfoEvidenceAnalysis::getVepSubFields, mismatch, 'vep' header sub field count: {}, parsed :{} sub fields"
                           , vep_header_ptr->subFieldHeaders().size(), vep_sub_fields.size());
      return std::nullopt;

    }

    if (std::string(vep_sub_fields[0]) == variant.alternate().getSequenceAsString()) {

      variant_vep_vector.emplace_back(std::move(vep_field));
      vep_sub_fields_vector.emplace_back(std::move(vep_sub_fields));

    }

  }

  // No point in returning an empty vep object.
  if (variant_vep_vector.empty()) {

    return std::nullopt;

  }

  // Return the object.
  return std::make_unique<const kgl::VEPSubFieldEvidence>(vep_header_ptr, std::move(variant_vep_vector), std::move(vep_sub_fields_vector));

}


// A temporary function to explore the vep field values.
// Returns all the distinct field values of a vep sub field from a population.
void kgl::InfoEvidenceAnalysis::vepSubFieldValues( std::string vep_sub_field,
                                                   const std::shared_ptr<const UnphasedPopulation>& population) {


  struct SubFieldValues{

    std::set<std::string> field_value_set_;
    std::string vep_sub_field_;

    bool getSubFieldValues(const std::shared_ptr<const Variant>& variant_ptr) {

      std::optional<std::unique_ptr<const kgl::VEPSubFieldEvidence>> vep_fields_opt = getVepSubFields(*variant_ptr);

      if (vep_fields_opt) {

        std::optional<size_t> vep_index_opt = vep_fields_opt.value()->vepHeader()->getSubFieldIndex(vep_sub_field_);

        if (not vep_index_opt) {

          ExecEnv::log().error("InfoEvidenceAnalysis::vepSubFieldValues, sub field: {} not found", vep_sub_field_);
          return false;

        }

        for (auto const& sub_fields : vep_fields_opt.value()->vepSubFields()) {

          field_value_set_.emplace(std::string(sub_fields[vep_index_opt.value()]));

        }

      }

      return true;

    }

  };

  SubFieldValues sub_field_values;
  sub_field_values.vep_sub_field_ = vep_sub_field;

  population->processAll(sub_field_values, &SubFieldValues::getSubFieldValues);

  for (auto const& field_value : sub_field_values.field_value_set_) {

    ExecEnv::log().info("vep field: {} has value: {}", vep_sub_field, field_value);

  }

}





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold parsed "vep" sub fields.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



