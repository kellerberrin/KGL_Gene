//
// Created by kellerberrin on 22/04/23.
//

#include "kgl_variant_factory_vcf_evidence_vep.h"



namespace kgl = kellerberrin::genome;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// vep sub-fields

const std::vector<std::string_view> kgl::VEPSubFieldEvidence::vepSubFields(const std::string& vep_field) {

  return Utility::viewTokenizer(vep_field, VEPSubFieldHeader::VEP_DELIMITER_CHAR);

}



// Creates a vep subfield object, if defined by the variant.
std::optional<std::unique_ptr<const kgl::VEPSubFieldEvidence>>
kgl::InfoEvidenceAnalysis::getVepSubFields(const Variant& variant) {

  // Check that the vep field exists for this variant.
  std::optional<const kgl::InfoSubscribedField> vep_field_opt = getSubscribedField( variant, VEPSubFieldHeader::VEP_FIELD_ID);

  if (not vep_field_opt) {

    return std::nullopt;

  }

  auto& vep_field_obj = vep_field_opt.value();

  // Get the vep header object.
  std::optional<std::shared_ptr<const VEPSubFieldHeader>> vep_header_opt = vep_field_obj.vepSubFieldHeader();

  if (not vep_header_opt) {

    ExecEnv::log().error("InfoEvidenceAnalysis::getVepSubFields, found 'vep' field but no matching 'vep' header defined");
    return std::nullopt;

  }

  std::shared_ptr<const VEPSubFieldHeader> vep_header_ptr = vep_header_opt.value();

  if (vep_header_ptr->subFieldHeaders().empty()) {

    ExecEnv::log().error("InfoEvidenceAnalysis::getVepSubFields, found empty 'vep' header");
    return std::nullopt;

  }

  if (not variant.evidence().infoData()) {

    ExecEnv::log().error("InfoEvidenceAnalysis::getVepSubFields, INFO data not defined for variant");
    return std::nullopt;

  }

  // Access the data block and retrieve the vep field vector.
  const DataMemoryBlock &data_block = *(variant.evidence().infoData().value());

  const std::vector<std::string> vep_field_vector = varianttoStrings(vep_field_obj.getData(data_block));

  // If the vector is empty then return nullopt.
  if (vep_field_vector.empty()) {

    ExecEnv::log().warn( "InfoEvidenceAnalysis::getVepSubFields; mismatch, Variant: {} has 'vep' header but no vep data",
                         variant.output(',', VariantOutputIndex::START_0_BASED, true));
    return std::nullopt;

  }

  std::vector<std::string> checked_field_vector;
  for (const auto& const_vep_field : vep_field_vector) {

    ++vep_records_;

    auto vep_field = const_vep_field;

    const std::vector<std::string_view>& sub_fields = VEPSubFieldEvidence::vepSubFields(vep_field);

    if (sub_fields.size() !=  vep_header_ptr->subFieldHeaders().size()) {

      // **** Gnomad 3 VEP bug workaround ****
      // This is added to handle the Gnomad 3 VEP bug where sub-sub-fields in the 'LoF_xxx' fields are delimited with ','
      // and are thus parsed as separate VEP records.

      // Issue a warning that we have encountered the Gnomad 3 problem.
      ++vep_size_errors_;
      static bool gnomad3_vep_warning{false};
      if (not gnomad3_vep_warning) {

        ExecEnv::log().warn("InfoEvidenceAnalysis::getVepSubFields; Gnomad 3 VEP bug, VEP sub-field count: {} not equal to VEP header size: {}",
                            sub_fields.size(),  vep_header_ptr->subFieldHeaders().size());
        ExecEnv::log().warn("InfoEvidenceAnalysis::getVepSubFields; Unparsed Vep Field: {}", vep_field);
        for (const auto& display_vep_field : vep_field_vector) {

          const std::vector<std::string_view>& display_sub_fields = VEPSubFieldEvidence::vepSubFields(display_vep_field);
          ExecEnv::log().warn("InfoEvidenceAnalysis::getVepSubFields; Vep Field: {}, Size: {}", display_vep_field, display_sub_fields.size());

        }
        gnomad3_vep_warning = true;

      }

    } else {

      checked_field_vector.push_back(std::move(vep_field));

    }

  }

  if (checked_field_vector.empty()) {

    return std::nullopt;

  }

  return std::make_unique<const kgl::VEPSubFieldEvidence>(vep_header_ptr, std::move(checked_field_vector));

}


kgl::VepValueMap kgl::InfoEvidenceAnalysis::getVepValues(const Variant& variant, std::vector<std::string> vep_field_list) {

  return getVepData(variant, getVepIndexes(variant, vep_field_list));

}

// Empty field list returns all available fields and associated offsets.
kgl::VepIndexVector kgl::InfoEvidenceAnalysis::getVepIndexes(const Variant& variant, const std::vector<std::string>& vep_field_list) {

  std::optional<std::unique_ptr<const VEPSubFieldEvidence>> vep_fields_opt = InfoEvidenceAnalysis::getVepSubFields(variant);

  if (not vep_fields_opt) {

    return VepIndexVector{};

  }

  auto const &vep_fields = *vep_fields_opt.value();
  auto const &vep_header = *vep_fields.vepHeader();

  VepIndexVector field_index;
  if (not vep_field_list.empty()) {

    for (auto const &field_header : vep_field_list) {

      std::optional<size_t> vep_index_opt = vep_header.getSubFieldIndex(field_header);

      if (not vep_index_opt) {

        ExecEnv::log().error("InfoEvidenceAnalysis::getVepValues, sub field: {} not found", field_header);
        return VepIndexVector{};

      }

      size_t vep_index = vep_index_opt.value();
      field_index.emplace_back(field_header, vep_index);

    }

  } else {

    size_t index{0};
    for (auto const &vep_field_header : vep_header.subFieldHeaders()) {

      field_index.emplace_back(vep_field_header, index);
      ++index;

    }

  }

  return field_index;

}

kgl::VepValueMap kgl::InfoEvidenceAnalysis::getVepData(const Variant& variant, const VepIndexVector& vep_field_list) {

  if (vep_field_list.empty()) {

    return VepValueMap{};

  }

  std::optional<std::unique_ptr<const VEPSubFieldEvidence>> vep_fields_opt = InfoEvidenceAnalysis::getVepSubFields(variant);

  if (not vep_fields_opt) {

    return VepValueMap{};

  }

  auto const &vep_fields = *vep_fields_opt.value();
  auto const &vep_header = *vep_fields.vepHeader();

  VepValueMap vep_value_map;
  for (auto const &vep_field : vep_fields.vepFields()) {

    const std::vector<std::string_view> sub_field_vector = VEPSubFieldEvidence::vepSubFields(vep_field);

    if (sub_field_vector.size() != vep_header.subFieldHeaders().size()) {

      ExecEnv::log().error("InfoEvidenceAnalysis::getVepData; Vep Header size: {} not equal vep sub field size: {}, vep field: {}",
                           vep_header.subFieldHeaders().size(), sub_field_vector.size(), vep_field);
      for (auto const &sub_field : sub_field_vector) {

        ExecEnv::log().error("InfoEvidenceAnalysis::getVepData; sub field: {}", std::string(sub_field));

      }

      return VepValueMap{};

    }

    std::map<std::string, std::string> sub_field_map;
    for (auto const&[field_ident, index] : vep_field_list) {

      sub_field_map[field_ident] = sub_field_vector[index];

    }

    vep_value_map.push_back(std::move(sub_field_map));

  }

  return vep_value_map;

}


// A temporary function to explore the vep field values.
// Returns all the distinct field values of a vep sub field from a population.
void kgl::InfoEvidenceAnalysis::vepSubFieldValues( std::string vep_sub_field,
                                                   const std::shared_ptr<const PopulationDB>& population) {


  VepSubFieldValues sub_field_values(vep_sub_field);

  population->processAll(sub_field_values, &VepSubFieldValues::getSubFieldValues);

  for (auto const& field_value : sub_field_values.getMap()) {

    ExecEnv::log().info("vep field: {} has value: {} count: {}", vep_sub_field, field_value.first, field_value.second);

  }

}


// A temporary function to explore the vep field values.
// Returns all the distinct field values of a vep sub field from a population.
void kgl::InfoEvidenceAnalysis::vepSubFieldCount( std::string vep_sub_field,
                                                  const std::shared_ptr<const PopulationDB>& population) {


  VepSubFieldValues sub_field_values(vep_sub_field);

  population->processAll(sub_field_values, &VepSubFieldValues::getSubFieldValues);

  ExecEnv::log().info("vep field: {} has count: {} discrete values", vep_sub_field, sub_field_values.getMap().size());

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::VepSubFieldValues::getSubFieldValues(const std::shared_ptr<const Variant>& variant_ptr) {

  std::optional<std::unique_ptr<const VEPSubFieldEvidence>> vep_fields_opt = InfoEvidenceAnalysis::getVepSubFields(*variant_ptr);

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
      if (field_value.empty()) {

        continue;

      }

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


// For a population
bool kgl::VepSubFieldValues::getPopulationValues(const std::shared_ptr<const PopulationDB>& population_ptr) {

  field_value_map_.clear();
  return population_ptr->processAll(*this, &VepSubFieldValues::getSubFieldValues);

}

// For a genome
bool kgl::VepSubFieldValues::getGenomeValues(const std::shared_ptr<const GenomeDB>& genome_ptr) {

  field_value_map_.clear();
  return genome_ptr->processAll(*this, &VepSubFieldValues::getSubFieldValues);

}

// For a contig.
bool kgl::VepSubFieldValues::getContigValues(const std::shared_ptr<const ContigDB>& contig_ptr) {

  field_value_map_.clear();
  return contig_ptr->processAll(*this, &VepSubFieldValues::getSubFieldValues);

}
