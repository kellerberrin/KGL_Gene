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


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Creates the static data and header indexes, and dynamically creates the InfoDataBlock object per VCF record.


void kgl::ManageInfoData::setupStaticStorage(InfoEvidenceHeader& evidence_header) {

  size_t field_index = 0;
  for (auto& subscribed_item : evidence_header.info_subscribed_map_) {

    InfoEvidenceIntern internal_type = subscribed_item.second.evidenceType().InternalInfoType();
    subscribed_item.second.fieldAddress(static_storage_.staticIncrementAndAllocate(internal_type));
    subscribed_item.second.fieldIndexId(field_index);
    ++field_index;

  }

}


// Create the storage to be used for the parsed info record.
std::unique_ptr<kgl::InfoDataBlock> kgl::ManageInfoData::setupDynamicStorage( const VCFInfoParser& info_parser,
                                                                              std::shared_ptr<const InfoEvidenceHeader> header_ptr) const {

// Parse the VCF info line.

  InfoDataUsageCount dynamic_info_storage = staticStorage();

  for (auto const &[ident, subscribed_info] : header_ptr->getMap()) {

    std::optional<InfoParserToken> token = info_parser.getToken(ident);

    // No additional storage required if token not available.
    if (token) {

      InfoEvidenceIntern internal_type = subscribed_info.evidenceType().InternalInfoType();
      if (not dynamic_info_storage.dynamicIncrementAndAllocate(internal_type, token.value())) {

        ExecEnv::log().warn("InfoDataUsageCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}, field ID:{}, Number:{}, Type:{}"
        , std::string(token.value().first), token.value().second, subscribed_info.infoVCF().ID,
                            subscribed_info.infoVCF().number, subscribed_info.infoVCF().type);

      }

    }

  }

  std::unique_ptr<InfoDataBlock> data_block_ptr(std::make_unique<InfoDataBlock>(header_ptr));

  data_block_ptr->allocateMemory(dynamic_info_storage);

  return std::move(data_block_ptr);

}



std::unique_ptr<kgl::InfoDataBlock> kgl::ManageInfoData::setupAndLoad( const VCFInfoParser& info_parser,
                                                                       std::shared_ptr<const InfoEvidenceHeader> self_ptr) const {

  std::unique_ptr<kgl::InfoDataBlock> data_block_ptr = setupDynamicStorage(info_parser, self_ptr);

  InfoDataUsageCount dynamic_accounting = staticStorage();  // Start with the static storage (indexes) already counted
  InfoDataUsageCount static_accounting; // Make sure this matches the object in staticStorage()
  for (auto& [ident, subscribed_item] : self_ptr->getMap()) {

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

  return manage_info_data_.setupAndLoad(info_parser, info_evidence_header_);

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
  manage_info_data_.setupStaticStorage(*info_evidence_header_);
  // Print out.
  std::string all_available = info_evidence_header_->getMap().size() == all_available_map_.size() ? "(all available)" : "";
  ExecEnv::log().info("Subscribed to {} {} VCF Info fields", info_evidence_header_->getMap().size(), all_available);

}


