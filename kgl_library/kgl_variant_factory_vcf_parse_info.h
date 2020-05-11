//
// Created by kellerberrin on 10/5/20.
//

#ifndef KGL_VARIANT_FACTORY_GRCH_INFO_H
#define KGL_VARIANT_FACTORY_GRCH_INFO_H

#include "kel_exec_env.h"
#include "kgl_variant.h"
#include "kgl_variant_evidence.h"

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



// Defines the automatically generated INFO field definitions for Gnomad VCF files.
struct VCFInfoDescription {

  bool process;
  std::string description;
  std::string type;
  std::string number;

};

using VCFInfoProcessMap = std::map<std::string, VCFInfoDescription>;


class GnomadEvidence : public VariantEvidence {

public:

  GnomadEvidence(size_t vcf_record_count) : VariantEvidence(vcf_record_count) {}
  ~GnomadEvidence() = default;

private:

  std::array<unsigned long, 20> test_{0};


};

// The efficient info parser and all the info field definitions for
// the GRCh38.r3.0 series of VCF files.
class GnomadInfo_3_0 {

public:

  explicit GnomadInfo_3_0(std::string&& vcf_info) : info_parser_(std::move(vcf_info)) {

    loadActive();

  }
  ~GnomadInfo_3_0() = default;

  [[nodiscard]] const VCFInfoParser& infoParser() const { return info_parser_; }

private:

  VCFInfoParser info_parser_;
  VCFInfoProcessMap active_map_;
  const static VCFInfoProcessMap gnomad_3_0_map_;

  // load the active INFO fields into the active_map.
  void loadActive();

};


// The efficient info parser and all the info field definitions for
// the GRCh38.r2.1.1 series of VCF files.
class GnomadInfo_2_1_1 {

public:

  explicit GnomadInfo_2_1_1(std::string&& vcf_info) : info_parser_(std::move(vcf_info)) {}
  ~GnomadInfo_2_1_1() = default;

  [[nodiscard]] const VCFInfoParser& infoParser() const { return info_parser_; }

private:

  VCFInfoParser info_parser_;
  VCFInfoProcessMap active_map_;
  const static VCFInfoProcessMap gnomad_2_1_1_map;


};



} // namespace




#endif //KGL_KGL_VARIANT_FACTORY_GRCH_INFO_H
