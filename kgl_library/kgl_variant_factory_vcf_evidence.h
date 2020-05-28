//
// Created by kellerberrin on 16/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_H


#include "kel_exec_env.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_vcf_evidence_data_mem.h"

#include <string>
#include <array>
#include <variant>


namespace kellerberrin::genome {   //  organization level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
using InfoDataVariant = std::variant<bool, std::vector<int64_t>, std::vector<double>, std::vector<std::string>, std::monostate>;
class ManageInfoData;
class InfoEvidenceHeader; // forward.

// Holds the indexing data to access the InfoDataBlock object and return VCF Info Data.
class InfoSubscribedField {

public:

//  friend ManageInfoData;

  InfoSubscribedField(VCFInfoRecord vcfInfoRecord,
                      std::shared_ptr<const InfoEvidenceHeader> info_evidence_header,
                      ManageInfoData& manage_info_data)
  : vcfInfoRecord_(std::move(vcfInfoRecord)),
    type_(InfoTypeLookup::evidenceType(vcfInfoRecord_)),
    info_evidence_header_(std::move(info_evidence_header)),
    m_data_handle_(requestResourceHandle(manage_info_data)) {}

  InfoSubscribedField(const InfoSubscribedField &) = default;
  ~InfoSubscribedField() = default;

  // Returns the original VCF header record for the VCF Info field.
  [[nodiscard]] const VCFInfoRecord &infoVCF() const { return vcfInfoRecord_; }
  // Returns one of the following
  // The enums are integer valued and can be used as an index for the std::variant return type of InfoDataVariant.
  // InfoEvidenceExtern::Boolean,
  // InfoEvidenceExtern::IntegerArray,
  // InfoEvidenceExtern::FloatArray,
  // InfoEvidenceExtern::StringArray,
  [[nodiscard]] InfoEvidenceExtern dataType() const { return evidenceType().ExternalInfoType(); }
  // Index into the variant type.
  [[nodiscard]] size_t dataIndex() const { return static_cast<size_t>(dataType()); }
  // The array size that will be returned. Zero if missing data.
  [[nodiscard]]  size_t dataSize(const InfoDataBlock& info_data_block) const {
    return info_data_block.getDataSize(fieldAddress(), fieldIndexId(), evidenceType().InternalInfoType()); }
  // Returns a std::variant containing the data which then can be indexed by the data type enum above.
  // Array is zero sized when there was no data item available.
  [[nodiscard]]  InfoDataVariant getData(const InfoDataBlock&) const;
  // The header object that can access all subscribed Info data items.
  std::shared_ptr<const InfoEvidenceHeader> getDataHeader() const { return info_evidence_header_; }

  // Need to limit visibility to these members.
  void fieldAddress(size_t field_address) { field_address_ = field_address; }
  void fieldIndexId(size_t field_index) { field_index_ = field_index; }
  [[nodiscard]] InfoEvidenceType evidenceType() const { return type_; }
  [[nodiscard]] size_t fieldAddress() const { return field_address_; }
  [[nodiscard]] size_t fieldIndexId() const { return field_index_; }

  const std::optional<InfoResourceHandle>& getDataHandle() const { return m_data_handle_; }

private:


  const VCFInfoRecord vcfInfoRecord_;  // The original VCF Header Record.
  const InfoEvidenceType type_;  // THe inferred subscriber type, external type and internal type.
  std::shared_ptr<const InfoEvidenceHeader> info_evidence_header_; // Ensure the index knows which header it belongs to.
  size_t field_address_{0};  // Permanent index into the InfoDataBlock object.
  size_t field_index_{0};  // Unique id integer.
  const std::optional<InfoResourceHandle> m_data_handle_;

  std::optional<InfoResourceHandle> requestResourceHandle(ManageInfoData& manage_info_data);

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An indexed map of InfoSubscribedFields. There is only one of these held by all variants with INFO evidence fields.

using InfoSubscribedMap = std::map<std::string, InfoSubscribedField>;
class EvidenceFactory;
class InfoEvidenceHeader {

public:

//  friend ManageInfoData;
  friend EvidenceFactory;

  explicit InfoEvidenceHeader() {}
  InfoEvidenceHeader(const InfoEvidenceHeader &) = delete;

  ~InfoEvidenceHeader() = default;

  [[nodiscard]] std::optional<const InfoSubscribedField> getSubscribedField(const std::string &field_id) const;
  [[nodiscard]] InfoSubscribedMap &getMap()  { return info_subscribed_map_; }
  [[nodiscard]] const InfoSubscribedMap &getConstMap() const { return info_subscribed_map_; }

private:

  InfoSubscribedMap info_subscribed_map_;

  // Note that this routine takes a shared_ptr to itself obtained from the info data factory. This is passed onto subscribed field objects.
  [[nodiscard]] bool setupEvidenceHeader( const VCFInfoRecord& vcf_info_record,
                                          std::shared_ptr<const InfoEvidenceHeader> self_ptr,
                                          ManageInfoData& manage_info_data);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Creates the static data and header indexes, and dynamically creates the InfoDataBlock object per VCF record.
// If we ever decide to change the way Info data is parsed and managed then we can swap out this object for another in the factory below.

class ManageInfoData {

public:

  ManageInfoData() {}
  ~ManageInfoData() = default;

  // Call this immediately after the evidence header has been created.
  void setupStaticStorage(InfoEvidenceHeader& evidence_header);
  // After the evidence header has been created, set variable ids, and addresses.
  [[nodiscard]] std::unique_ptr<InfoDataBlock> setupAndLoad( const VCFInfoParser& info_parser,
                                                             std::shared_ptr<const InfoEvidenceHeader> evidence_ptr) const;


  [[nodiscard]] InfoMemoryResource& resourceAllocator() { return resource_allocator_; }


private:

  [[nodiscard]] const InfoDataUsageCount& staticStorage() const { return static_storage_; }
  [[nodiscard]] std::unique_ptr<InfoDataBlock> setupDynamicStorage( const VCFInfoParser& info_parser,
                                                                    std::shared_ptr<const InfoEvidenceHeader> evidence_ptr) const;

  InfoDataUsageCount static_storage_;

  // Initialized on creation of the EvidenceHeader.
  InfoMemoryResource resource_allocator_;

  // Copy the resource allocator and resolve all dynamic resource requests by examining the parsed data.
  // Return the dynamically resolved resource object.
  [[nodiscard]] InfoMemoryResource resolveResources(const VCFInfoParser& info_parser, const InfoEvidenceHeader& evidence_header) const;

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
  // All subscribed Info fields.
  std::shared_ptr<const InfoEvidenceHeader> getInfoHeader() const { return info_evidence_header_; }


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

  // Remove after testing.
  void debugData(const VCFInfoParser& info_parser, std::unique_ptr<InfoDataBlock>& info_data_ptr) const;

};



} // namespace





#endif //KGL_VARIANT_FACTORY_VCF_EVIDENCE_H
