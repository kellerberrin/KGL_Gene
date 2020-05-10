//
// Created by kellerberrin on 10/5/20.
//

#ifndef KGL_VARIANT_FACTORY_GRCH_INFO_H
#define KGL_VARIANT_FACTORY_GRCH_INFO_H

#include "kgl_variant_factory_vcf_parse_impl.h"

#include <string>
#include <string_view>


namespace kellerberrin::genome {   //  organization level namespace



// Defines the automatically INFO field definitions for Gnomad VCF files.
struct VCFInfoDescription {

  bool process;
  std::string description;
  std::string type;
  std::string number;

};

using VCFInfoProcessMap = std::map<std::string, VCFInfoDescription>;

// The efficient info parser and all the info field definitions for
// the GRCh38.r3.0 series of VCF files.
class GnomadInfo_3_0 {

public:

  explicit GnomadInfo_3_0(std::string&& vcf_info) : info_parser_(std::move(vcf_info)) {}
  ~GnomadInfo_3_0() = default;

  [[nodiscard]] const VCFInfoParser& infoParser() const { return info_parser_; }

private:

  VCFInfoParser info_parser_;
  const static VCFInfoProcessMap gnomad_3_0_map;


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
  const static VCFInfoProcessMap gnomad_2_1_1_map;


};






} // namespace




#endif //KGL_KGL_VARIANT_FACTORY_GRCH_INFO_H
