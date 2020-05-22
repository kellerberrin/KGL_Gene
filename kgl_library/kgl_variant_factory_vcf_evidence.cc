//
// Created by kellerberrin on 16/5/20.
//

#include "kgl_variant_factory_vcf_evidence.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_vcf_parse_info.h"

namespace kgl = kellerberrin::genome;




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper class for INFO types.


// Returns true if a number > 1. Used in the type definition lambdas.
bool kgl::InfoTypeLookup::isVectorType(const std::string& type) {

  try {

    return std::stol(type) > 1;

  }
  catch(...) {

    return false;

  }

}


kgl::InfoEvidenceType kgl::InfoTypeLookup::evidenceType(const VCFInfoRecord &vcf_info_item) {

  for (auto const& evidence_type : type_definitions_) {

    if (evidence_type.second(vcf_info_item.type, vcf_info_item.number)) {

      return evidence_type.first;

    }

  }

  ExecEnv::log().warn("InfoTypeLookup::evidenceType, Info ID: {} Unable to find data type combination for Number: {}, Type: {}",
                      vcf_info_item.ID, vcf_info_item.number, vcf_info_item.type);

  return {InfoEvidenceSubscriber::NotImplemented, InfoEvidenceExtern::NotImplemented, InfoEvidenceIntern::NotImplemented };

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

bool kgl::InfoEvidenceType::fixedDataType() const {

  switch (InternalInfoType()) {

    case InfoEvidenceIntern::intern_char: // Size is fixed (boolean variables) and known at Info data subscription time.
    case InfoEvidenceIntern::intern_integer: // Size is fixed and known at Info data subscription time.
    case InfoEvidenceIntern::intern_float: // Size is fixed and known at Info data subscription time.
    case InfoEvidenceIntern::intern_integer_array: // Size is fixed and known at Info data subscription time.
    case InfoEvidenceIntern::intern_float_array: // Size is fixed and known at Info data subscription time.
    case InfoEvidenceIntern::NotImplemented:  // Trivially fixed.
      return true;

    case InfoEvidenceIntern::intern_string: // Size varies between records.
    case InfoEvidenceIntern::intern_string_array:  // Size varies between records.
    case InfoEvidenceIntern::intern_unity_integer_array:   // Size varies between records.
    case InfoEvidenceIntern::intern_unity_float_array:   // Size varies between records.
    case InfoEvidenceIntern::intern_unity_string_array:    // Size varies between records.
      return false;

    default:
      return false;
  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// data_item_count is the number of data items, for a float array of size 3 this is 1 for 1 array. For a string array of size 3 it is 3.
// data_item_size is the number of underlying data items. For float array of size 3 this is 3.
// For a string array this will be the total number of characters in all strings. For a string array of size 3 data_item_count=3, but data_item_size is size of
// all three strings added together. The total characters required for all 3 strings.
kgl::ItemOffset kgl::DataInfoTypeCount::staticIncrementAndAllocate(InfoEvidenceIntern internal_type) {

  ItemOffset item_offset;
  switch (internal_type) {

    case InfoEvidenceIntern::intern_char:  {

      item_offset.offset = char_count_;
      item_offset.is_array = false;
      ++char_count_;

    }// Size is fixed (boolean variables) and known at Info data subscription time.
    return item_offset;

    case InfoEvidenceIntern::intern_unity_integer_array:
    case InfoEvidenceIntern::intern_integer: {

      item_offset.offset = integer_count_;
      item_offset.is_array = false;
      ++integer_count_;

    }// Size is fixed and known at Info data subscription time.
    return item_offset;

    case InfoEvidenceIntern::intern_unity_float_array:
    case InfoEvidenceIntern::intern_float: {

      item_offset.offset = float_count_;
      item_offset.is_array = false;
      ++float_count_;  // should be 1

    }// Size is fixed and known at Info data subscription time.
    return item_offset;

    case InfoEvidenceIntern::intern_string: {

      item_offset.offset = string_count_;
      item_offset.is_array = false;
      ++string_count_;   // allocate a std::string_view

    } // Size varies between records.
    return item_offset;

    case InfoEvidenceIntern::intern_integer_array: {

      item_offset.offset = array_count_;
      item_offset.is_array = true;
      ++array_count_;

    } // Size is fixed and known at Info data subscription time.
    return item_offset;

    case InfoEvidenceIntern::intern_float_array: {

      item_offset.offset = array_count_;
      item_offset.is_array = true;
      ++array_count_;

    }
    return item_offset;


    case InfoEvidenceIntern::intern_string_array: {

      item_offset.offset = array_count_;
      item_offset.is_array = true;
      ++array_count_;

    } // Size varies between records.
    return item_offset;

    case InfoEvidenceIntern::intern_unity_string_array: {

      item_offset.offset = array_count_;
      item_offset.is_array = true;
      ++array_count_;

    }   // Size varies between records.
    return item_offset;

    case InfoEvidenceIntern::NotImplemented:  // Trivially fixed.
    default:
      return item_offset;

  }

}


bool kgl::DataInfoTypeCount::dynamicIncrementAndAllocate(const InfoSubscribedField& subscribed_field, const InfoParserToken& token) {

  InfoEvidenceIntern internal_type = subscribed_field.evidenceType().InternalInfoType();

  switch (internal_type) {

    // Data is pre-allocated for the fixed fields.
      case InfoEvidenceIntern::intern_char:
        if (token.second != 0) {

          ExecEnv::log().warn("DataInfoTypeCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}, field ID:{}, Number:{}, Type:{}"
          , std::string(token.first), token.second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);
          return false;
        }
        break;

      case InfoEvidenceIntern::intern_integer:
      case InfoEvidenceIntern::intern_float:
        if (token.second != 1) {

          ExecEnv::log().warn("DataInfoTypeCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}, field ID:{}, Number:{}, Type:{}"
          , std::string(token.first), token.second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);
          return false;
        }
        break;

      case InfoEvidenceIntern::intern_unity_integer_array:
        if (token.second > 1) {

          ++unity_array_count_;
          integer_count_ += token.second;

        }
        break;

      case InfoEvidenceIntern::intern_unity_float_array:
        if (token.second > 1) {

          ++unity_array_count_;
          float_count_ += token.second;

        }
        break;

      case InfoEvidenceIntern::intern_string: {

        char_count_ += token.first.size(); // total char size of all strings.

      }
      break;

      case InfoEvidenceIntern::intern_integer_array: {

        integer_count_ += token.second; // size of the array

      } // Size is fixed and known at Info data subscription time.
      break;

      case InfoEvidenceIntern::intern_float_array: {

        float_count_ += token.second; // size of the array

      }
      break;

      case InfoEvidenceIntern::intern_unity_string_array:
      case InfoEvidenceIntern::intern_string_array: {

        char_count_ += token.first.size() - (token.second - 1); // total char size of all strings, less the delimiter chars.
        string_count_ += token.second;   // number strings, allocate a vector of std::string_views

      } // Size varies between records.
      break;

     case InfoEvidenceIntern::NotImplemented:  // Trivially fixed.
     default:
      break;

  }

  return true;

}


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

  for (auto& subscribed_item : info_subscribed_map_) {

    InfoEvidenceIntern internal_type = subscribed_item.second.evidenceType().InternalInfoType();
    subscribed_item.second.dataOffset(static_storage_.staticIncrementAndAllocate(internal_type));

  }

}

// Create the storage to be used for the parsed info record.
std::unique_ptr<kgl::InfoDataBlock> kgl::InfoEvidenceHeader::setupDynamicStorage( const VCFInfoParser& info_parser,
                                                                                  std::shared_ptr<const InfoEvidenceHeader> self_ptr) const {

// Parse the VCF info line.

  DataInfoTypeCount dynamic_info_storage = staticStorage();

  for (auto const &subscribed_info : getMap()) {

    std::optional<InfoParserToken> token = info_parser.getToken(subscribed_info.first);

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

  return info_evidence_header_->setupDynamicStorage(info_parser, info_evidence_header_);

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
