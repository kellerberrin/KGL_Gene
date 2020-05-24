//
// Created by kellerberrin on 16/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_H


#include "kel_exec_env.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_vcf_evidence_memory.h"

#include <string>
#include <array>


namespace kellerberrin::genome {   //  organization level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

class ManageInfoData;
class InfoEvidenceHeader; // forward.
// Holds the indexing data to access the InfoDataBlock object and return VCF Info Data.
class InfoSubscribedField {

public:

  friend ManageInfoData;

  InfoSubscribedField(VCFInfoRecord vcfInfoRecord, std::shared_ptr<const InfoEvidenceHeader> info_evidence_header)
  : vcfInfoRecord_(std::move(vcfInfoRecord)),
    type_(InfoTypeLookup::evidenceType(vcfInfoRecord_)),
    info_evidence_header_(std::move(info_evidence_header)) {}
  InfoSubscribedField(const InfoSubscribedField &) = default;
  ~InfoSubscribedField() = default;

  [[nodiscard]] const VCFInfoRecord &infoRecord() const { return vcfInfoRecord_; }
  [[nodiscard]]  InfoEvidenceExtern dataType() const { return evidenceType().ExternalInfoType(); }


private:

  void fieldAddress(size_t field_address) { field_address_ = field_address; }
  void fieldIndexId(size_t field_index) { field_index_ = field_index; }
  [[nodiscard]] size_t fieldAddress() const { return field_address_; }
  [[nodiscard]] size_t fieldIndexId() const { return field_index_; }
  [[nodiscard]] InfoEvidenceType evidenceType() const { return type_; }

  const VCFInfoRecord vcfInfoRecord_;  // The original VCF Header Record.
  const InfoEvidenceType type_;  // THe inferred subscriber type, external type and internal type.
  std::shared_ptr<const InfoEvidenceHeader> info_evidence_header_; // Ensure the index knows which header it belongs to.
  size_t field_address_{0};  // Permanent index into the InfoDataBlock object.
  size_t field_index_{0};  // Unique id integer.

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An indexed map of InfoSubscribedFields. There is only one of these held by all variants with INFO evidence fields.

using InfoSubscribedMap = std::map<std::string, InfoSubscribedField>;
class EvidenceFactory;
class InfoEvidenceHeader {

public:

  friend ManageInfoData;
  friend EvidenceFactory;

  explicit InfoEvidenceHeader() {}
  InfoEvidenceHeader(const InfoEvidenceHeader &) = delete;

  ~InfoEvidenceHeader() = default;

  [[nodiscard]] std::optional<const InfoSubscribedField> getSubscribedField(const std::string &field_id) const;
  [[nodiscard]] const InfoSubscribedMap &getMap() const { return info_subscribed_map_; }


private:

  InfoSubscribedMap info_subscribed_map_;

  // Note that this routine takes a shared_ptr to itself obtained from the info data factory. This is passed onto subscribed field objects.
  [[nodiscard]] bool setupEvidenceHeader(const VCFInfoRecord& vcf_info_record, std::shared_ptr<const InfoEvidenceHeader> self_ptr);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Creates the static data and header indexes, and dynamically creates the InfoDataBlock object per VCF record.
// If we ever decide to change the way Info data is parsed and managed then we can swap out this object for another in the factory below.

class ManageInfoData {

public:

  ManageInfoData() = default;
  ~ManageInfoData() = default;

  // Call this immediately after the evidence header has been created.
  void setupStaticStorage(InfoEvidenceHeader& evidence_header);
  // After the evidence header has been created, set variable ids, and addresses.
  [[nodiscard]] std::unique_ptr<InfoDataBlock> setupAndLoad( const VCFInfoParser& info_parser, std::shared_ptr<const InfoEvidenceHeader> self_ptr) const;

private:

  [[nodiscard]] const InfoDataUsageCount& staticStorage() const { return static_storage_; }
  [[nodiscard]] std::unique_ptr<InfoDataBlock> setupDynamicStorage(const VCFInfoParser& info_parser, std::shared_ptr<const InfoEvidenceHeader> self_ptr) const;

  InfoDataUsageCount static_storage_;

  // Static storage and field offsets. Indexes the fields.
  // Fixed storage that does not change betweeen info records.

};




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The evidence factory creates a common evidence lookup object for all variants (to optimize memory usage).
// The evidence factory also creates an evidence object for each variant.

// createVariantEvidence() either returns a single data block, InfoDataBlock for a single alternate allele or
// a MultipleAlleleDataBlock data block for a multiple alternative allele VCF record.
using InfoDataEvidence = std::optional<std::shared_ptr<InfoDataBlock>>;
class EvidenceFactory {

public:

  explicit EvidenceFactory(const EvidenceInfoSet& evidence_set)  // An ordered set of Info field identifiers.
  : evidence_map_(evidence_set),
    info_evidence_header_(std::make_shared<InfoEvidenceHeader>()) {}
  ~EvidenceFactory() = default;

  // A map of all available Info fields parsed from the VCF header
  // This also initializes the InfoEvidenceHeader object.
  void availableInfoFields(const VCFInfoRecordMap& vcf_info_map);
  [[nodiscard]] const VCFInfoRecordMap& availableInfoFields() const { return all_available_map_; }
  // For each input VCFRecord info text field (std::moved for efficiency), create a parsed data object.
  [[nodiscard]] InfoDataEvidence createVariantEvidence(std::string&& info);


private:

  // The evidence fields specified in the runtime XML file.
  const EvidenceInfoSet evidence_map_;
  // All available info fields.
  VCFInfoRecordMap all_available_map_;
  // The Info header block, contains all the subscribed field information.
  std::shared_ptr<InfoEvidenceHeader> info_evidence_header_;
  // Manage the definition and creation of Info Data objects.
  ManageInfoData manage_info_data_;

  // If the user specifies just specifies "None" (case insensitive) then no Info fields will be subscribed.
  constexpr static const char *NO_FIELD_SUBSCRIBED_ = "NONE";

};



} // namespace





#endif //KGL_VARIANT_FACTORY_VCF_EVIDENCE_H
