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

    std::vector<std::string> bin_strings = Utility::charTokenizer(bin_data.front(), BIN_DELIMITER_);

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

    const auto& info_data = *(variant.evidence().infoData().value());
    return info_data.evidenceHeader()->getSubscribedField(field_ident);

  }

  return std::nullopt;

}


std::optional<kgl::InfoDataVariant> kgl::InfoEvidenceAnalysis::getInfoData( const Variant& variant,
                                                                            const std::string& field_ident) {

  std::optional<const kgl::InfoSubscribedField> field_opt = getSubscribedField(variant, field_ident);

  if (field_opt) {

    auto info_data_ptr = variant.evidence().infoData();
    if (info_data_ptr) {

      const DataMemoryBlock &data_block = *(info_data_ptr.value());

      InfoDataVariant variant_data = field_opt->getData(data_block);

      return variant_data;

    }

  }

  return std::nullopt;

}

