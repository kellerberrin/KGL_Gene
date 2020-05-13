//
// Created by kellerberrin on 10/5/20.
//

#ifndef KGL_VARIANT_FACTORY_PARSE_INFO_H
#define KGL_VARIANT_FACTORY_PARSE_INFO_H

#include "kel_exec_env.h"
#include "kgl_variant.h"
#include "kgl_variant_evidence.h"
#include "kgl_variant_factory_vcf_parse_header.h"

#include <string>
#include <string_view>
#include <array>


namespace kellerberrin::genome {   //  organization level namespace




/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is an efficient single pass parser.
// In general std::string_view is work of the devil and a seg fault waiting to happen.
// But if the underlying string has the same guaranteed lifetime as the associated std::string_view(s) then a seg fault
// may not be inevitable.
using InfoParserMap = std::map<std::string_view, std::string_view>;
class VCFInfoParser {

public:

  // std::move the info string into this object for efficiency.
  explicit VCFInfoParser(std::string&& info) : info_(std::move(info)), info_view_(info_) {

    if (not parseInfo()) {

      ExecEnv::log().error("VCFInfoParser::VCFInfoParser, Problem parsing info field");

    }

  }
  ~VCFInfoParser() = default;

  [[nodiscard]] const InfoParserMap& getMap() const { return parsed_token_map_; }
  [[nodiscard]] const std::string& info() const { return info_; }
  [[nodiscard]] std::optional<std::string> getInfoField(const std::string& key) const;

private:

  std::string info_;
  std::string_view info_view_;
  InfoParserMap parsed_token_map_;

  constexpr static const char INFO_FIELD_DELIMITER_{';'};
  constexpr static const char INFO_VALUE_DELIMITER_{'='};

  [[nodiscard]] bool parseInfo();

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The evidence factory creates a common evidence lookup object for all variants (to optimize memory usage).
// The evidence factory also creates an evidence object for each variant.


class EvidenceFactory {

public:

  explicit EvidenceFactory(const EvidenceInfoSet& evidence_map) : evidence_map_(evidence_map) {}
  ~EvidenceFactory() = default;

  void availableInfoFields(const VCFInfoRecordMap& vcf_info_map);
  [[nodiscard]] std::unique_ptr<VariantEvidence> createVariantEvidence(size_t record_count, std::string&& info);

private:

  // The evidence fields specified in the runtime XML file.
  const EvidenceInfoSet evidence_map_;
  // available info fields parsed from the VCF header.
  VCFInfoRecordMap active_info_map_;

};



} // namespace




#endif //KGL_KGL_VARIANT_FACTORY_GRCH_INFO_H
