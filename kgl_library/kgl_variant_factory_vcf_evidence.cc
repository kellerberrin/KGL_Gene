//
// Created by kellerberrin on 16/5/20.
//

#include "kgl_variant_factory_vcf_evidence.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_vcf_evidence_data_blk.h"

#include <variant>

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//


//  Returns the field size, A size of zero is a missing field.
kgl::InfoDataVariant kgl::InfoSubscribedField::getData(const InfoDataBlock& info_data_block) const {

  switch (evidenceType().InternalInfoType()) {

    // Data is pre-allocated for the fixed fields.
    case InfoEvidenceIntern::intern_char:
      return info_data_block.getBoolean(fieldAddress());

    case InfoEvidenceIntern::intern_integer: {

      std::vector<int64_t> int_array;
      if (info_data_block.getInteger(fieldAddress()) != InfoDataBlock::MISSING_VALUE_INTEGER_) {

        int_array.emplace_back(static_cast<int64_t>(info_data_block.getInteger(fieldAddress())));

      }
      return int_array;

    }

    case InfoEvidenceIntern::intern_float: {

      std::vector<double> float_array;

      if (info_data_block.getFloat(fieldAddress()) != InfoDataBlock::MISSING_VALUE_FLOAT_) {

        float_array.emplace_back(static_cast<double>(info_data_block.getFloat(fieldAddress())));

      }
      return float_array;

    }

    case InfoEvidenceIntern::intern_unity_integer_array:
      return info_data_block.getUnityIntegerArray(fieldAddress(), fieldIndexId());

    case InfoEvidenceIntern::intern_unity_float_array:
      return info_data_block.getUnityFloatArray(fieldAddress(), fieldIndexId());

    case InfoEvidenceIntern::intern_string: {

      std::vector<std::string> string_array;
      string_array.emplace_back(info_data_block.getString(fieldAddress()));
      return string_array;

    }

    case InfoEvidenceIntern::intern_integer_array:
      return info_data_block.getIntegerArray(fieldAddress(), fieldIndexId());

    case InfoEvidenceIntern::intern_float_array:
      return info_data_block.getFloatArray(fieldAddress(), fieldIndexId());

    case InfoEvidenceIntern::intern_unity_string_array:
    case InfoEvidenceIntern::intern_string_array:
      return info_data_block.getStringArray(fieldAddress(), fieldIndexId());

    case InfoEvidenceIntern::NotImplemented:  // unknown internal type.
    default:
      ExecEnv::log().error( "InfoDataBlock::getDataSize, Internal data type unknown, cannot get size");
      return std::monostate();
  }

}

//  Returns the data, A size of zero is a missing field.
kgl::InfoDataVariant kgl::InfoSubscribedField::getNewData(const DataMemoryBlock& memory_block) const {

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
      ExecEnv::log().error( "InfoSubscribedField::getNewData, Internal data type unknown");
      return std::monostate();

  }

}


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
// An indexed map of InfoSubscribedField. There is only one of these held by all variants with INFO evidence fields.

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


std::optional<const kgl::InfoSubscribedField> kgl::InfoEvidenceHeader::getSubscribedField(const std::string& field_id) const {

  auto result = info_subscribed_map_.find(field_id);

  if (result == info_subscribed_map_.end()) {

    ExecEnv::log().warn("InfoEvidenceHeader::getSubscribedField, Info field ID : {} is not a subscribed field", field_id);
    return std::nullopt;

  }


  return result->second;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Creates the static data and header indexes, and dynamically creates the InfoDataBlock object per VCF record.


void kgl::ManageInfoData::setupStaticStorage(InfoEvidenceHeader& evidence_header) {

  size_t field_index = 0;
  for (auto& subscribed_item : evidence_header.getMap()) {

    InfoEvidenceIntern internal_type = subscribed_item.second.evidenceType().InternalInfoType();
    subscribed_item.second.fieldAddress(static_storage_.staticIncrementAndAllocate(internal_type));
    subscribed_item.second.fieldIndexId(field_index);
    ++field_index;

  }

}


// Create the storage to be used for the parsed info record.
std::unique_ptr<kgl::InfoDataBlock> kgl::ManageInfoData::setupDynamicStorage( const VCFInfoParser& info_parser,
                                                                              std::shared_ptr<const InfoEvidenceHeader> evidence_ptr) const {

// Parse the VCF info line.

  InfoDataUsageCount dynamic_info_storage = staticStorage();

  for (auto const &[ident, subscribed_info] : evidence_ptr->getConstMap()) {

    std::optional<InfoParserToken> token = info_parser.getToken(ident);

    // No additional storage required if token not available.
    if (token) {

      InfoEvidenceIntern internal_type = subscribed_info.evidenceType().InternalInfoType();
      if (not dynamic_info_storage.dynamicIncrementAndAllocate(subscribed_info.getDataHandle().initialDataSize(), internal_type, token.value())) {

        ExecEnv::log().warn("InfoDataUsageCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}, field ID:{}, Number:{}, Type:{}"
        , std::string(token.value().first), token.value().second, subscribed_info.infoVCF().ID,
                            subscribed_info.infoVCF().number, subscribed_info.infoVCF().type);

      }

    }

  }

  std::unique_ptr<InfoDataBlock> data_block_ptr(std::make_unique<InfoDataBlock>(evidence_ptr));

  data_block_ptr->allocateMemory(dynamic_info_storage);

  return std::move(data_block_ptr);

}

std::unique_ptr<const kgl::DataMemoryBlock> kgl::ManageInfoData::createMemoryBlock( const VCFInfoParser& info_parser,
                                                                                    std::shared_ptr<const InfoEvidenceHeader> evidence_ptr) const {

  InfoMemoryResource resolved_resource = resolveResources(info_parser, *evidence_ptr);

  return std::make_unique<const DataMemoryBlock>(evidence_ptr, resolved_resource, info_parser);

}

std::unique_ptr<kgl::InfoDataBlock> kgl::ManageInfoData::setupAndLoad( const VCFInfoParser& info_parser,
                                                                       std::shared_ptr<const InfoEvidenceHeader> evidence_ptr) const {



  std::unique_ptr<kgl::InfoDataBlock> data_block_ptr = setupDynamicStorage(info_parser, evidence_ptr);


  InfoDataUsageCount dynamic_accounting = staticStorage();  // Start with the static storage (indexes) already counted
  InfoDataUsageCount static_accounting; // Make sure this matches the object in staticStorage()
  for (auto& [ident, subscribed_item] : evidence_ptr->getConstMap()) {

    const std::optional<InfoParserToken>& parser_token = info_parser.getToken(ident);

    if (not data_block_ptr->indexAndVerify(subscribed_item.fieldAddress(),
                                           subscribed_item.fieldIndexId(),
                                           subscribed_item.evidenceType().InternalInfoType(),
                                           parser_token,
                                           dynamic_accounting,
                                           static_accounting)) {

      // Any error from this function is a memory block error, so report the error and terminate.
      std::string token_value = parser_token ? std::string(parser_token.value().first) : "MISSING_VALUE";
      size_t token_size = parser_token ? parser_token.value().second : 0;
      ExecEnv::log().error("ManageInfoData::setupAndLoad, Problem loading data for VCF Info field value: {}, size:{}", token_value, token_size);
      ExecEnv::log().error("ManageInfoData::setupAndLoad, VCF Info field ID: {}, Number: {}, Type: {}, Description: {}",
                           subscribed_item.infoVCF().ID, subscribed_item.infoVCF().number,
                           subscribed_item.infoVCF().type,
                           subscribed_item.infoVCF().description);
      ExecEnv::log().critical("ManageInfoData::setupAndLoad, Variant VCF Info field memory block set up encountered a serious error and cannot continue ...");

    }

  }

  // Final checks for memory usage.
  // This checks that the raw data allocated exactly matches the data utilized by all the Info fields.
  if (not (dynamic_accounting == data_block_ptr->getRawMemoryUsage())) {

    ExecEnv::log().error("ManageInfoData::setupAndLoad, The Dynamic Accounting Object Not Equal to the Total Allocated Data size");
    ExecEnv::log().critical("ManageInfoData::setupAndLoad, Variant VCF Info field memory block set up encountered a serious error and cannot continue ...");

  }

  // This checks that the data allocated for field indexing exactly matches the indexes used by all Info fields.
  if (not (static_accounting == staticStorage())) {

    ExecEnv::log().error("ManageInfoData::setupAndLoad, The Static Accounting Object Not Equal to the Static Data size");
    ExecEnv::log().critical("ManageInfoData::setupAndLoad, Variant VCF Info field memory block set up encountered a serious error and cannot continue ...");

  }

  // Data block containing all subscribed Info fields.
  return std::move(data_block_ptr);

}


kgl::InfoMemoryResource kgl::ManageInfoData::resolveResources( const VCFInfoParser& info_parser,
                                                               const InfoEvidenceHeader& evidence_header) const {

  // Copy the static resource allocator.
  InfoMemoryResource resolved_resources = resource_allocator_;

  for (auto const &[ident, subscribed_info] : evidence_header.getConstMap()) {

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


kgl::InfoDataEvidence kgl::EvidenceFactory::createVariantEvidence(std::string&& info) {

  // If no Info fields have been subscribed, then just return std::nullopt
  if (info_evidence_header_->getMap().empty()) {

    return std::nullopt;

  }

  // Parse the VCF info line.
  VCFInfoParser info_parser(std::move(info));

  std::unique_ptr<InfoDataBlock> info_data_ptr = manage_info_data_.setupAndLoad(info_parser, info_evidence_header_);

  std::unique_ptr<const DataMemoryBlock> mem_blk_ptr = manage_info_data_.createMemoryBlock(info_parser, info_evidence_header_);

  // Debug code
  debugData(info_parser, info_data_ptr, mem_blk_ptr);

  return std::move(mem_blk_ptr);

}


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

  } else { // subscribe to the specified fields.

    for (auto const &ident : evidence_map_) {

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

  // Setup the static storage allocation and field offsets.
  manage_info_data_.setupStaticStorage(*info_evidence_header_);
  // Print out.
  std::string all_available = info_evidence_header_->getMap().size() == all_available_map_.size() ? "(all available)" : "";
  ExecEnv::log().info("Subscribed to {} {} VCF Info fields", info_evidence_header_->getMap().size(), all_available);

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Debug stuff. To be removed.

void kgl::EvidenceFactory::debugData( const VCFInfoParser& info_parser,
                                      std::unique_ptr<InfoDataBlock>& info_data_ptr,
                                      std::unique_ptr<const DataMemoryBlock>& mem_blk_ptr) const {

  for (auto const& [ident, info_item] : getInfoHeader()->getConstMap()) {

    InfoDataVariant item_data = info_item.getData(*info_data_ptr);
    if (item_data.index() != info_item.dataIndex()) {

      ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Info Field: {},  std::variant index: {} not equal to return index: {}"
      , ident, item_data.index(), info_item.dataIndex());

    }

  }

  for (auto const& [ident, field_item] : info_evidence_header_->getMap()) {

    size_t data_size = field_item.dataSize(*info_data_ptr);

    auto token = info_parser.getToken(ident);
////////////////////////////////////////////////
    if (token) {

      if (field_item.dataType() != InfoEvidenceExtern::Boolean and token.value().second != data_size) {

        if (field_item.dataType() != InfoEvidenceExtern::FloatArray and std::string(token.value().first) != ".") {

          ExecEnv::log().error("EvidenceFactory::createVariantEvidence, Parser Count not equal InfoDataBlock Field: {}, Number: {}, Type: {}, Info Size: {}, Parser Size: {} Value: {}",
                               field_item.infoVCF().ID, field_item.infoVCF().number, field_item.infoVCF().type, data_size, token.value().second, std::string(token.value().first));

        }

        switch (field_item.dataType()) {

          case InfoEvidenceExtern::Boolean:
            if (not std::get<bool>(field_item.getData(*info_data_ptr))) {

              ExecEnv::log().error("EvidenceFactory::createVariantEvidence, Boolean should be true: {}, Number: {}, Type: {}, Info Size: {}, Parser Size: {} Value: {}",
                                   field_item.infoVCF().ID, field_item.infoVCF().number, field_item.infoVCF().type, data_size, token.value().second, std::string(token.value().first));
            }
            break;

          case InfoEvidenceExtern::IntegerArray: {

            std::vector<int64_t> int_vector = std::get<std::vector<int64_t>>(field_item.getData(*info_data_ptr));
            InfoParserIntegerArray integer_array = VCFInfoParser::getInfoIntegerArray(std::string(token.value().first));

            if (int_vector.size() == integer_array.size()) {

              for (size_t index = 0; index < integer_array.size(); ++index) {

                if (int_vector[index] != integer_array[index]) {

                  ExecEnv::log().error(
                  "EvidenceFactory::createVariantEvidence, Info Integer: {} not equal Token integer: {}, Id: {}, Number: {}, Type: {}, Info Size: {}, Parser Size: {} Value: {}",
                  int_vector[index], integer_array[index], field_item.infoVCF().ID, field_item.infoVCF().number,
                  field_item.infoVCF().type, data_size, token.value().second, std::string(token.value().first));


                }

              }

            } else if (not( int_vector.empty() and std::string(token.value().first) == ".")) {

              ExecEnv::log().error(
              "EvidenceFactory::createVariantEvidence, Info Integer Size: {} not equal Token Integer Size: {}, Id: {}, Number: {}, Type: {}, Info Size: {}, Parser Size: {} Value: {}",
              int_vector.size(), integer_array.size(), field_item.infoVCF().ID, field_item.infoVCF().number,
              field_item.infoVCF().type, data_size, token.value().second, std::string(token.value().first));

            }

          }
            break;

          case InfoEvidenceExtern::StringArray: {

            std::vector<std::string> string_vector = std::get<std::vector<std::string>>(field_item.getData(*info_data_ptr));
            std::vector<std::string_view> string_view_vector = VCFInfoParser::getInfoStringArray(token.value().first);

            if (string_vector.size() == string_view_vector.size()) {

              for (size_t index = 0; index < string_vector.size(); ++index) {

                if (string_vector[index] != std::string(string_view_vector[index])) {

                  ExecEnv::log().error(
                  "EvidenceFactory::createVariantEvidence, Info String: {} not equal Token string: {}, Id: {}, Number: {}, Type: {}, Info Size: {}, Parser Size: {} Value: {}",
                  string_vector[index], std::string(string_view_vector[index]), field_item.infoVCF().ID,
                  field_item.infoVCF().number, field_item.infoVCF().type, data_size, token.value().second,
                  std::string(token.value().first));


                }

              }

            } else  {

              ExecEnv::log().error(
              "EvidenceFactory::createVariantEvidence, Info String Size: {} not equal Token String Size: {}, Id: {}, Number: {}, Type: {}, Info Size: {}, Parser Size: {} Value: {}",
              string_vector.size(), string_view_vector.size(), field_item.infoVCF().ID,
              field_item.infoVCF().number, field_item.infoVCF().type, data_size, token.value().second,
              std::string(token.value().first));

            }

          }
            break;

          case InfoEvidenceExtern::FloatArray: {

            std::vector<double> float_vector = std::get<std::vector<double>>(field_item.getData(*info_data_ptr));
            InfoParserFloatArray float_array =  VCFInfoParser::getInfoFloatArray(std::string(token.value().first));

            if (float_vector.size() == float_array.size()) {

              for (size_t index = 0; index < float_array.size(); ++index) {

                if (float_vector[index] != float_array[index]) {

                  ExecEnv::log().error(
                  "EvidenceFactory::createVariantEvidence, Info Float: {} not equal Token Float: {}, Id: {}, Number: {}, Type: {}, Info Size: {}, Parser Size: {} Value: {}",
                  float_vector[index], float_array[index], field_item.infoVCF().ID, field_item.infoVCF().number,
                  field_item.infoVCF().type, data_size, token.value().second, std::string(token.value().first));


                }

              }

            } else if (not(float_vector.empty() and std::string(token.value().first) == ".")) {

              ExecEnv::log().error(
              "EvidenceFactory::createVariantEvidence, Info Float Size: {} not equal Token Float Size: {}, Id: {}, Number: {}, Type: {}, Info Size: {}, Parser Size: {} Value: {}",
              float_vector.size(), float_array.size(), field_item.infoVCF().ID, field_item.infoVCF().number,
              field_item.infoVCF().type, data_size, token.value().second, std::string(token.value().first));


            }

          }
            break;

          case InfoEvidenceExtern::NotImplemented:
            ExecEnv::log().error("EvidenceFactory::createVariantEvidence, Type not Implemented for, Field: {}, Number: {}, Type: {}, Info Size: {}",
                                 field_item.infoVCF().ID, field_item.infoVCF().number, field_item.infoVCF().type, data_size);
            break;

        }

      }

    } else if (data_size != 0) {

      if (field_item.dataType() != InfoEvidenceExtern::Boolean or std::get<bool>(field_item.getData(*info_data_ptr)))

        ExecEnv::log().error("EvidenceFactory::createVariantEvidence, Present in InfoData but absent from Parser, Field: {}, Number: {}, Type: {}, Info Size: {}",
                             field_item.infoVCF().ID, field_item.infoVCF().number, field_item.infoVCF().type, data_size);

    }
////////////////////////////////////////////////////////
  }


  for (auto const& [ident, info_item] : getInfoHeader()->getConstMap()) {


    InfoDataVariant item_data = info_item.getData(*info_data_ptr);
    InfoDataVariant new_item_data = info_item.getNewData(*mem_blk_ptr);
    if (item_data != new_item_data) {

      // check the variant index.

      if (item_data.index() != new_item_data.index()) {

        ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Info Field: {}, handle: {}  info data index: {} not equal to mem blk index: {}"
        , ident, info_item.getDataHandle().handleId(), item_data.index(), new_item_data.index());

      } else { // check values.

        switch(info_item.dataType()) {

          case InfoEvidenceExtern::Boolean:
            if (std::get<bool>(item_data) != std::get<bool>(new_item_data)) {

              ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Boolean difference Old: {} New: {}"
              , std::get<bool>(item_data) ? "True" : "False", std::get<bool>(new_item_data) ? "True" : "False");

            }
            break;

          case InfoEvidenceExtern::FloatArray:
            if (std::get<std::vector<double>>(item_data).size() != std::get<std::vector<double>>(new_item_data).size()) {

              ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Ident: {}, Float Difference Vector Size, Old: {} New: {}"
              , ident, std::get<std::vector<double>>(item_data).size(), std::get<std::vector<double>>(new_item_data).size());

              if (std::get<std::vector<double>>(item_data).size() == 1) {

                ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Ident: {}, Float Value[0], Old: {} "
                , ident, std::get<std::vector<double>>(item_data)[0]);

              }


            } else {

              size_t index = 0;
              for (auto float_val : std::get<std::vector<double>>(item_data)) {

                if (float_val != std::get<std::vector<double>>(new_item_data)[index]) {

                  auto token = info_parser.getToken(ident);
                  if (token) {

                    ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Float Vector element mismatch ident: {}, index: {}, Old: {} New: {}, Token: {}"
                    , ident, index, float_val, std::get<std::vector<double>>(new_item_data)[index], std::string(token.value().first));

                  } else {

                    ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Float Vector element mismatch ident: {}, index: {}, Old: {} New: {}, NoToken"
                    , ident, index, float_val, std::get<std::vector<double>>(new_item_data)[index]);

                  }

                }

                ++index;

              }

            }
            break;

          case InfoEvidenceExtern::IntegerArray:
            if (std::get<std::vector<int64_t>>(item_data).size() != std::get<std::vector<int64_t>>(new_item_data).size()) {

              ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Integer Difference Vector Size, Old: {} New: {}"
              , std::get<std::vector<int64_t>>(item_data).size(), std::get<std::vector<int64_t>>(new_item_data).size());

            } else {

              size_t index = 0;
              for (auto integer_val : std::get<std::vector<int64_t>>(item_data)) {

                if (integer_val != std::get<std::vector<int64_t>>(new_item_data)[index]) {

                  ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Integer Vector element mismatch index: {}, Old: {} New: {}"
                  , index, integer_val, std::get<std::vector<int64_t>>(new_item_data)[index]);

                }

                ++index;

              }

            }
            break;

          case InfoEvidenceExtern::StringArray:
            if (std::get<std::vector<std::string>>(item_data).size() != std::get<std::vector<std::string>>(new_item_data).size()) {

              ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, String Difference Vector Size, Old: {} New: {}"
              , std::get<std::vector<std::string>>(item_data).size(), std::get<std::vector<std::string>>(new_item_data).size());

            } else {

              size_t index = 0;
              for (auto string_val : std::get<std::vector<std::string>>(item_data)) {

                if (string_val != std::get<std::vector<std::string>>(new_item_data)[index]) {

                  ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, String Vector element mismatch index: {}, Old: {} New: {}"
                  , index, string_val, std::get<std::vector<std::string>>(new_item_data)[index]);

                }

                ++index;

              }

            }
            break;

          case InfoEvidenceExtern::NotImplemented:
            ExecEnv::log().error("InfoEvidenceHeader::debugReturnAll, Data type is 'Not Implemented'");
            break;

        }

      }

    }

  }

}


