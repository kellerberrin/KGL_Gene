//
// Created by kellerberrin on 16/5/20.
//

#include "kgl_variant_factory_vcf_evidence.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_vcf_evidence_data_blk.h"

#include "kel_utility.h"

#include <variant>

namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Extract VEP subfields.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The VEP field is an Gnomad specific field that has (generally 68) numerous sub-fields.
// For convenience these can be accessed by name.
// Non-Gnomad Info data will simply return a std::nullopt.

bool kgl::VEPSubFieldHeader::parseHeader(const std::string& description) {

  size_t header_offset = description.find(HEADER_SEARCH_STR_);

  if (header_offset == std::string::npos) {

    ExecEnv::log().error("VEPSubFieldHeader::parseHeader, could not parse vep header string: {}", description);
    return false;

  }

  header_offset += std::string(HEADER_SEARCH_STR_).size();

  std::string unparsed_header = description.substr(header_offset);

  sub_fields_headers_ = Utility::charTokenizer(unparsed_header, VEP_DELIMITER_CHAR);

  size_t index{0};
  for (auto const& sub_field : sub_fields_headers_) {

    auto result = index_map_.try_emplace(sub_field, index);

    if (not result.second) {

      ExecEnv::log().error("VEPSubFieldHeader::parseHeader, duplicate sub field header text: {}", sub_field);
      return false;

    }

    ++index;

  }

  return true;

}


std::optional<size_t> kgl::VEPSubFieldHeader::getSubFieldIndex(const std::string& sub_field) const {

  auto result = index_map_.find(sub_field);

  if (result == index_map_.end()) {

    return std::nullopt;

  }

  return result->second;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Returns the data as a variant of vectors and a boolean, an empty vector is a missing field.
//
kgl::InfoDataVariant kgl::InfoSubscribedField::getData(const DataMemoryBlock& memory_block) const {

  // Check that this is the correct index for the data block.
  // Should never happen (but best to check).
  if (info_evidence_header_ != memory_block.evidenceHeader()) {

    ExecEnv::log().error( "InfoSubscribedField::getData, Incorrect Subscribed Field Index used Access Info Data Block");
    return std::monostate(); // Inaccessible variant data, will trigger downstream errors.

  }

  switch(getDataHandle().resourceType()) {

    case DataResourceType::Boolean:
      return memory_block.getBoolean(getDataHandle());

    case DataResourceType::Integer:
      return memory_block.getInteger(getDataHandle());

    case DataResourceType::Float:
      return memory_block.getFloat(getDataHandle());

    case DataResourceType::String:
      return memory_block.getString(getDataHandle());

    default:
      ExecEnv::log().error( "InfoSubscribedField::getData, Internal data type unknown");
      return std::monostate(); // Inaccessible variant data, will trigger upstream errors.

  }

}

// Manage internally (compressed) stored data types and map onto external (std::variant) data types.
kgl::InfoResourceHandle kgl::InfoSubscribedField::requestResourceHandle(ManageInfoData& manage_info_data) {

  InfoMemoryResource& resource_allocator = manage_info_data.resourceAllocator();

  switch (evidenceType().InternalInfoType()) {

    // Data is pre-allocated for the fixed fields.
    case InfoEvidenceIntern::intern_char:
      return resource_allocator.resourceRequest(DataResourceType::Boolean, DataDynamicType::FixedData, 1);

    case InfoEvidenceIntern::intern_integer:
      return resource_allocator.resourceRequest(DataResourceType::Integer, DataDynamicType::FixedData, 1);

    case InfoEvidenceIntern::intern_float:
      return resource_allocator.resourceRequest(DataResourceType::Float, DataDynamicType::FixedData, 1);

    case InfoEvidenceIntern::intern_unity_integer_array:
      return resource_allocator.resourceRequest(DataResourceType::Integer, DataDynamicType::FixedDynamic, 1);

    case InfoEvidenceIntern::intern_unity_float_array:
      return resource_allocator.resourceRequest(DataResourceType::Float, DataDynamicType::FixedDynamic, 1);

    case InfoEvidenceIntern::intern_string:
      return resource_allocator.resourceRequest(DataResourceType::String, DataDynamicType::FixedData, 1);

    case InfoEvidenceIntern::intern_integer_array:
      return resource_allocator.resourceRequest(DataResourceType::Integer, DataDynamicType::DynamicData, 0);

    case InfoEvidenceIntern::intern_float_array:
      return resource_allocator.resourceRequest(DataResourceType::Float, DataDynamicType::DynamicData, 0);

    case InfoEvidenceIntern::intern_unity_string_array:
    case InfoEvidenceIntern::intern_string_array:
      return resource_allocator.resourceRequest(DataResourceType::String, DataDynamicType::DynamicData, 0);

    case InfoEvidenceIntern::NotImplemented:  // unknown internal type.
    default:
      ExecEnv::log().error( "InfoSubscribedField::requestResourceHandle, Internal data type unknown, cannot obtain data handle");
      // Just return a boolean descriptor.
      return resource_allocator.resourceRequest(DataResourceType::Boolean, DataDynamicType::FixedData, 1);
  }

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An indexed map of InfoSubscribedFields. There is only one of these held by all variants with INFO evidence fields
// parsed from a particular VCF file. However, this fact is opaque to the library user, as all access is through
// the getSubscribedField() function below where fields are requested by name.

bool kgl::InfoEvidenceHeader::setupEvidenceHeader( const VCFInfoRecord& vcf_info_record,
                                                   std::shared_ptr<const InfoEvidenceHeader> self_ptr,
                                                   ManageInfoData& manage_info_data) {

  InfoSubscribedField subscribed_field(vcf_info_record, self_ptr, manage_info_data);
  auto result = info_subscribed_map_.try_emplace(vcf_info_record.ID, subscribed_field);

  if (not result.second) {

    ExecEnv::log().warn("InfoEvidenceHeader::setupEvidenceHeader. (Duplicate) unable to subscribe to duplicate info field: {}", vcf_info_record.ID);
    return false;

  }

  return true;

}

// All field access is through this function.
std::optional<const kgl::InfoSubscribedField> kgl::InfoEvidenceHeader::getSubscribedField(const std::string& field_id) const {

  auto result = info_subscribed_map_.find(field_id);

  if (result == info_subscribed_map_.end()) {

    return std::nullopt;

  }

  auto const [field_name, field_info] = *result;

  return field_info;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Creates the static data and header indexes, and dynamically creates the InfoDataBlock object per VCF record.

std::shared_ptr<const kgl::DataMemoryBlock> kgl::ManageInfoData::createMemoryBlock( const VCFInfoParser& info_parser,
                                                                                    std::shared_ptr<const InfoEvidenceHeader> evidence_ptr) const {

  InfoMemoryResource resolved_resource = resolveResources(info_parser, *evidence_ptr);

  return std::make_shared<const DataMemoryBlock>(evidence_ptr, resolved_resource, info_parser);

}



kgl::InfoMemoryResource kgl::ManageInfoData::resolveResources( const VCFInfoParser& info_parser,
                                                               const InfoEvidenceHeader& evidence_header) const {

  // Copy the static resource allocator.
  InfoMemoryResource resolved_resources = resource_allocator_;

  // For all subscribed Info fields.
  for (auto const &[ident, subscribed_info] : evidence_header.getConstMap()) {

    // The pre-parsed info field. This is used to determine storage requirements.
    std::optional<InfoParserToken> token = info_parser.getToken(ident);

    // No additional storage required if token not available.
    if (token) {

        // Resolve the runtime memory allocation by examining the data.
      if (not resolved_resources.resolveAllocation(subscribed_info.getDataHandle(), token.value())) {

        ExecEnv::log().warn("ManageInfoData::resolveResources, Bad size (expected 1) Token: {} size: {}, field ID:{}, Number:{}, Type:{}"
        , std::string(token.value().first), token.value().second, subscribed_info.infoVCF().ID,
                            subscribed_info.infoVCF().number, subscribed_info.infoVCF().type);

      }

    }

  }

  return resolved_resources;

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The evidence factory creates a common evidence lookup object for all variants (to optimize memory usage).
// The evidence factory also creates an evidence object for each variant (data only).

// Create an indexed data object from the VCF info field, std::moved in as a string.
std::shared_ptr<const kgl::DataMemoryBlock> kgl::EvidenceFactory::createVariantEvidence(std::string&& info) {

  // If no Info fields have been subscribed, then just return std::nullopt
  if (info_evidence_header_->getMap().empty()) {

    return nullptr;

  }

  // Parse the VCF info line.
  VCFInfoParser info_parser(std::move(info));

  // Use the parsed data to create a compact memory block with a copy of the Info data.
  std::shared_ptr<const DataMemoryBlock> mem_blk_ptr = manage_info_data_.createMemoryBlock(info_parser, info_evidence_header_);

  return mem_blk_ptr;

}

// Subscribe to info fields.
void kgl::EvidenceFactory::availableInfoFields(const VCFInfoRecordMap& vcf_info_map) {


  all_available_map_ = vcf_info_map;

  // If no info fields are specified then all available fields are subscribed.
  if (evidence_map_.empty()) {

    ExecEnv::log().info("No specific Info fields specified, subscribing to all available fields");

    for (auto const& [ident, vcf_info_record] : all_available_map_) {

      if (not info_evidence_header_->setupEvidenceHeader(vcf_info_record, info_evidence_header_, manage_info_data_)) {

        ExecEnv::log().warn("EvidenceFactory::availableInfoFields, Cannot subscribe to Info field: {} description :{}", ident, vcf_info_record.description);

      }

    }

  } else { // Subscribe to a list of explicitly specified fields.

    for (auto const &ident : evidence_map_) {

      // If the subscribed field is "None", then no fields are subscribed.
      // Note from above, if there are no fields are explicitly subscribed, then they are all subscribed.
      if (Utility::toupper(ident) == NO_FIELD_SUBSCRIBED_) {

        continue;  // skip

      }

      auto result = all_available_map_.find(ident);

      if (result == all_available_map_.end()) {

        ExecEnv::log().warn("InfoEvidenceHeader::setupEvidenceHeader. Unable to subscribe to Info field: {}, field not found in list of available fields", ident);

      }
      else {

        if (not info_evidence_header_->setupEvidenceHeader(result->second, info_evidence_header_, manage_info_data_)) {

          ExecEnv::log().warn("EvidenceFactory::availableInfoFields, (Duplicate) unable to subscribe to duplicate info field: {}", ident);

        }

      }

    }

  }

  // Print out subscribed fields.
  std::string all_available = info_evidence_header_->getMap().size() == all_available_map_.size() ? "(all available)" : "";
  ExecEnv::log().info("Subscribed to {} {} VCF Info fields", info_evidence_header_->getMap().size(), all_available);

}

