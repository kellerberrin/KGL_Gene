//
// Created by kellerberrin on 5/6/20.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"

#include <vector>

namespace kgl = kellerberrin::genome;


const std::vector<std::string>& kgl::InfoEvidenceAnalysis::varianttoStrings(const InfoDataVariant& info_data) {

  static const std::vector<std::string> empty_vector;
  auto p_string_vector = std::get_if<std::vector<std::string>>(&info_data);
  if (p_string_vector != nullptr) {

    return *p_string_vector;

  } else {

    ExecEnv::log().error("InfoEvidenceAnalysis::varianttoStrings, Info data variant does not contain a string vector.");
    return empty_vector;

  }

}


const std::vector<double>& kgl::InfoEvidenceAnalysis::varianttoFloats(const InfoDataVariant& info_data) {

  static const std::vector<double> empty_vector;
  auto p_float_vector = std::get_if<std::vector<double>>(&info_data);
  if (p_float_vector != nullptr) {

    return *p_float_vector;

  } else {

    ExecEnv::log().error("InfoEvidenceAnalysis::varianttoFloats, Info data variant does not contain a float vector.");
    return empty_vector;

  }

}



const std::vector<int64_t>& kgl::InfoEvidenceAnalysis::varianttoIntegers(const InfoDataVariant& info_data) {

  static const std::vector<int64_t> empty_vector;
  auto p_integer_vector = std::get_if<std::vector<int64_t>>(&info_data);
  if (p_integer_vector != nullptr) {

    return *p_integer_vector;

  } else {

    ExecEnv::log().error("InfoEvidenceAnalysis::varianttoIntegers, Info data variant does not contain a integer vector.");
    return empty_vector;

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


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// vep sub-fields

const std::vector<std::string_view> kgl::VEPSubFieldEvidence::vepSubFields(const std::string& vep_field) {

  return Utility::view_tokenizer(vep_field, VEPSubFieldHeader::VEP_DELIMITER_CHAR);

}



// Creates a vep subfield object, if defined by the variant.
std::optional<std::unique_ptr<const kgl::VEPSubFieldEvidence>>
kgl::InfoEvidenceAnalysis::getVepSubFields(const Variant& variant) {

  static size_t correct_parse{0};
  // Check that the vep field exists for this variant.
  std::optional<const kgl::InfoSubscribedField> vep_field_opt = getSubscribedField( variant, VEPSubFieldHeader::VEP_FIELD_ID);

  if (not vep_field_opt) {

    return std::nullopt;

  }

  // Get the vep header object.
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

  // Access the data block and retrieve the vep field vector.
  const DataMemoryBlock &data_block = *variant.evidence().infoData().value();

  std::vector<std::string> vep_field_vector = varianttoStrings(vep_field_opt.value().getData(data_block));

  // If the vector is empty then return nullopt.
  if (vep_field_vector.empty()) {

    ExecEnv::log().warn( "InfoEvidenceAnalysis::getVepSubFields; mismatch, Variant: {} has 'vep' header but no vep data",
                          variant.output(',', VariantOutputIndex::START_0_BASED, true));
    return std::nullopt;

  }

  for (auto const& vep_field : vep_field_vector) {

    const std::vector<std::string_view>& sub_fields = VEPSubFieldEvidence::vepSubFields(vep_field);

    if (sub_fields.size() !=  vep_header_ptr->subFieldHeaders().size()) {

      ExecEnv::log().error("InfoEvidenceAnalysis::getVepSubFields; VEP sub-field count: {} not equal to VEP header size: {}, vep field: {}, variants correctly parsed: {}",
                           sub_fields.size(),  vep_header_ptr->subFieldHeaders().size(), vep_field, correct_parse);

      for (auto const& sub_field : sub_fields) {

        ExecEnv::log().error("InfoEvidenceAnalysis::vepSubFieldValues; sub field: {}", std::string(sub_field));

      }
      return std::nullopt;

    }

  }

  ++correct_parse;
  return std::make_unique<const kgl::VEPSubFieldEvidence>(vep_header_ptr, std::move(vep_field_vector));

}


// A temporary function to explore the vep field values.
// Returns all the distinct field values of a vep sub field from a population.
void kgl::InfoEvidenceAnalysis::vepSubFieldValues( std::string vep_sub_field,
                                                   const std::shared_ptr<const PopulationVariant>& population) {

  struct SubFieldValues{

    std::map<std::string, size_t> field_value_map_;
    std::string vep_sub_field_;

    bool getSubFieldValues(const std::shared_ptr<const Variant>& variant_ptr) {

      std::optional<std::unique_ptr<const kgl::VEPSubFieldEvidence>> vep_fields_opt = getVepSubFields(*variant_ptr);

      if (vep_fields_opt) {

        auto const &vep_fields = *vep_fields_opt.value();
        auto const &vep_header = *vep_fields.vepHeader();
        std::optional<size_t> vep_index_opt = vep_header.getSubFieldIndex(vep_sub_field_);

        if (not vep_index_opt) {

          ExecEnv::log().error("InfoEvidenceAnalysis::vepSubFieldValues, sub field: {} not found", vep_sub_field_);
          return false;

        }

        size_t vep_index = vep_index_opt.value();

        for (auto const&  vep_field : vep_fields.vepFields()) {

          const std::vector<std::string_view> sub_field_vector = VEPSubFieldEvidence::vepSubFields(vep_field);

            // Check that the header and the field vector are the same size.
          if (sub_field_vector.size() != vep_header.subFieldHeaders().size()) {

            ExecEnv::log().error("InfoEvidenceAnalysis::vepSubFieldValues; Vep Header size: {} not equal vep sub field size: {}, vep field: {}",
                                 vep_header.subFieldHeaders().size(), sub_field_vector.size(), vep_field);
            for (auto const& sub_field : sub_field_vector) {

              ExecEnv::log().error("InfoEvidenceAnalysis::vepSubFieldValues; sub field: {}", std::string(sub_field));

            }
            return false;
          }

          std::string field_value(sub_field_vector[vep_index]);
          auto result = field_value_map_.find(field_value);

          if (result != field_value_map_.end()) {

            ++result->second;

          } else {

            auto insert_result = field_value_map_.try_emplace(field_value, 1);
            if (not insert_result.second) {

              ExecEnv::log().error("InfoEvidenceAnalysis::vepSubFieldValues, could not insert sub_field value: {}", field_value);

            }

          }

        } // sub fields

      } // for all sub fields

      return true;

    } // member function.

  }; // struct.

  SubFieldValues sub_field_values;
  sub_field_values.vep_sub_field_ = vep_sub_field;

  population->processAll(sub_field_values, &SubFieldValues::getSubFieldValues);

  for (auto const& field_value : sub_field_values.field_value_map_) {

    ExecEnv::log().info("vep field: {} has value: {} count: {}", vep_sub_field, field_value.first, field_value.second);

  }

}


