//
// Created by kellerberrin on 16/5/20.
//

#include "kgl_variant_factory_vcf_evidence.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_vcf_parse_info.h"


namespace kgl = kellerberrin::genome;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An indexed map of InfoSubscribedField. There is only one of these held by all variants with INFO evidence fields.

bool kgl::InfoEvidenceHeader::setupEvidenceHeader(const VCFInfoRecord& vcf_info_record, std::shared_ptr<const InfoEvidenceHeader> self_ptr) {

  InfoSubscribedField subscribed_field(vcf_info_record, self_ptr);
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


void kgl::InfoEvidenceHeader::setupStaticStorage() {

  size_t field_index = 0;
  for (auto& subscribed_item : info_subscribed_map_) {

    InfoEvidenceIntern internal_type = subscribed_item.second.evidenceType().InternalInfoType();
    subscribed_item.second.dataOffset(static_storage_.staticIncrementAndAllocate(internal_type));
    subscribed_item.second.fieldIndex(field_index);
    ++field_index;

  }

}

// Create the storage to be used for the parsed info record.
std::unique_ptr<kgl::InfoDataBlock> kgl::InfoEvidenceHeader::setupDynamicStorage( const VCFInfoParser& info_parser,
                                                                                  std::shared_ptr<const InfoEvidenceHeader> self_ptr) const {

// Parse the VCF info line.

  DataInfoTypeCount dynamic_info_storage = staticStorage();

  for (auto const &subscribed_info : getMap()) {

    std::optional<InfoParserToken> token = info_parser.getToken(subscribed_info.first);

    // No additional storage required if token not available.
    if (token) {

      if (not dynamic_info_storage.dynamicIncrementAndAllocate(subscribed_info.second, token.value())) {

        ExecEnv::log().warn("EvidenceFactory::parseSubscribed, Vector value for assumed Scalar Field");

      }

    }

  }

  std::unique_ptr<InfoDataBlock> data_block_ptr(std::make_unique<InfoDataBlock>(self_ptr));

  data_block_ptr->allocateMemory(dynamic_info_storage);

  return std::move(data_block_ptr);

}


std::unique_ptr<kgl::InfoDataBlock> kgl::InfoEvidenceHeader::setupAndLoad( const VCFInfoParser& info_parser,
                                                                           std::shared_ptr<const InfoEvidenceHeader> self_ptr) const {

  std::unique_ptr<kgl::InfoDataBlock> data_block_ptr = setupDynamicStorage(info_parser, self_ptr);

  DataInfoTypeCount dynamic_accounting = staticStorage();  // start with the static storage already counted
  DataInfoTypeCount static_accounting; // make sure this matches the object in staticStorage()
  std::vector<ItemOffset> index_accounting;  // make sure this matches the permanent indexes in the InfoSubscribedField objects.
  for (auto& [ident, subscribed_item] : info_subscribed_map_) {

    std::optional<InfoParserToken> token = info_parser.getToken(ident);

    if (not data_block_ptr->indexAndVerify(subscribed_item, token, dynamic_accounting, static_accounting, index_accounting)) {

      std::string token_value = token ? std::string(token.value().first) : "MISSING_VALUE";
      ExecEnv::log().error("Problem loading data for field: {}, value: {}", token_value);

    }

  }

  // Finally, the usage accounting object should equal the original memory allocated.

//  if (dynamic_accounting != data_block_ptr->getTypeCount()) {
//
//
//  }

  return std::move(data_block_ptr);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// InfoDataBlock

kgl::DataInfoTypeCount kgl::InfoDataBlockNaive::dataPayload() {


//  static std::mutex func_lock;
//  std::scoped_lock<std::mutex> raii(func_lock);


  size_t float_count = 0;
  size_t integer_count = 0;
  size_t char_count = 0;
  size_t array_count = 0;
  size_t string_count = 0;

  char_count = bool_data_.size();

  integer_count = integer_data_.size();

  float_count = float_data_.size();

  size_t char_bytes = 0;
  for (auto const& string_field :  string_data_) {

    ++string_count;
    char_bytes += string_field.size();

  }

  for (auto const& float_array :  float_array_data_) {

    if (float_array.size() <= 1) {

      ++float_count;

    } else {

      float_count += float_array.size();
      ++array_count;

    }

  }

  for (auto const&integer_array :  integer_array_data_) {

    if (integer_array.size() <= 1) {

      ++float_count;

    } else {

      integer_count += integer_array.size();
      ++array_count;

    }

  }

  for (auto const& string_array :  string_array_data_) {

    if (string_array_data_.size() == 1) {

      ++string_count;
      char_bytes += string_array_data_.front().size();

    } else if (string_array_data_.empty() ) {

      ++string_count;

    } else {

      ++array_count;

      for (auto const& string_item : string_array) {

        ++string_count;
        char_bytes += string_item.size();

      }

    }

  }

  DataInfoTypeCount info_count;

  info_count.integerCount(integer_count);
  info_count.floatCount(float_count);
  info_count.charCount(char_count);
  info_count.arrayCount(array_count);
  info_count.stringCount(string_count);

  return info_count;

}

void kgl::InfoDataBlockNaive::clearAll() {

// The stored data.
  bool_data_.clear();
  float_data_.clear();
  integer_data_.clear();
  string_data_.clear();
  float_array_data_.clear();
  integer_array_data_.clear();
  string_array_data_.clear();

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The evidence factory creates a common evidence lookup object for all variants (to optimize memory usage).
// The evidence factory also creates an evidence object for each variant (data only).


kgl::InfoDataEvidence kgl::EvidenceFactory::createVariantEvidence(std::string&& info) {

  // If no Info fields have been subscribed, then just return std::nullopt
  if (info_evidence_header_->getMap().empty()) {

    return std::nullopt;

  }

//  std::unique_ptr<InfoDataBlockNaive> temp = parseSubscribed(std::move(info));

//  DataInfoTypeCount type_count = );

//  return std::move(temp);

//  std::unique_ptr<InfoDataBlock> data_block_ptr(std::make_unique<InfoDataBlock>(info_evidence_header_));

//  data_block_ptr->allocateMemory(type_count);

  return parseSubscribed_alt(std::move(info));

}


std::unique_ptr<kgl::InfoDataBlock> kgl::EvidenceFactory::parseSubscribed_alt(std::string&& info) {

  // Parse the VCF info line.
  VCFInfoParser info_parser(std::move(info));

  return info_evidence_header_->setupAndLoad(info_parser, info_evidence_header_);

}

std::unique_ptr<kgl::InfoDataBlockNaive> kgl::EvidenceFactory::parseSubscribed(std::string&& info) {

  // Parse the VCF info line.
  VCFInfoParser info_parser(std::move(info));

  std::unique_ptr info_data_ptr(std::make_unique<InfoDataBlockNaive>(info_evidence_header_));
  // Fill all the subscribed Info field values.
  for (auto const& subscribed_info : info_evidence_header_->getMap()) {

    switch(subscribed_info.second.evidenceType().ExternalInfoType()) {

      case InfoEvidenceExtern::Float: {
        info_data_ptr->float_data_.emplace_back(info_parser.getInfoFloat(subscribed_info.first));
      }
      break;

      case InfoEvidenceExtern::Integer: {

        info_data_ptr->integer_data_.emplace_back(info_parser.getInfoInteger(subscribed_info.first));

      }
      break;

      case InfoEvidenceExtern::String: {

        info_data_ptr->string_data_.emplace_back(info_parser.getInfoString(subscribed_info.first));

      }
      break;

      case InfoEvidenceExtern::Boolean: {

        info_data_ptr->bool_data_.emplace_back(info_parser.getInfoBoolean(subscribed_info.first));

      }
      break;

      case InfoEvidenceExtern::IntegerArray: {

        info_data_ptr->integer_array_data_.emplace_back(info_parser.getInfoIntegerArray(subscribed_info.first));

      }
      break;

      case InfoEvidenceExtern::FloatArray: {

        info_data_ptr->float_array_data_.emplace_back(info_parser.getInfoFloatArray(subscribed_info.first));

      }
      break;

      case InfoEvidenceExtern::StringArray: {

        info_data_ptr->string_array_data_.emplace_back(info_parser.getInfoStringArray(subscribed_info.first));

      }
      break;

      default:
        ExecEnv::log().error("EvidenceFactory::parseSubscribed, Unexpected INFO data type for field: {}", subscribed_info.first);
        break;

    }

  }

  return info_data_ptr;

}


void kgl::EvidenceFactory::availableInfoFields(const VCFInfoRecordMap& vcf_info_map) {


  all_available_map_ = vcf_info_map;

// If no info fields are specified then all available fields are subscribed.
  if (evidence_map_.empty()) {

    ExecEnv::log().info("No specific Info fields specified, subscribing to all available fields");

    for (auto const& [ident, vcf_info_record] : all_available_map_) {

      if (not info_evidence_header_->setupEvidenceHeader(vcf_info_record, info_evidence_header_)) {

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

        if (not info_evidence_header_->setupEvidenceHeader(result->second, info_evidence_header_)) {

          ExecEnv::log().warn("EvidenceFactory::availableInfoFields, (Duplicate) unable to subscribe to duplicate info field: {}", ident);

        }

      }

    }

  }

  // Setup the static storage allocation and field offsets.
  info_evidence_header_->setupStaticStorage();
  // Print out.
  std::string all_available = info_evidence_header_->getMap().size() == all_available_map_.size() ? "(all available)" : "";
  ExecEnv::log().info("Subscribed to {} {} VCF Info fields", info_evidence_header_->getMap().size(), all_available);

}
