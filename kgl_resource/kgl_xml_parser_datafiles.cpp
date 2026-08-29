//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_xml_parsers.h"
#include "kgl_xml_reader.h"
#include "kgl_xml_detail.h"

#include <string>
#include <vector>


namespace kgl = kellerberrin::genome;
namespace kglx = kellerberrin::genome::xml;

namespace {

  constexpr std::string_view DATA_FILE_LIST_{"dataFileList"};
  constexpr std::string_view DATA_FILE_IDENT_{"dataFileIdent"};
  constexpr std::string_view DATA_FILE_NAME_{"dataFileName"};
  constexpr std::string_view DATA_PARSER_TYPE_{"dataFileType"};
  constexpr std::string_view VCF_DATA_FILE_TYPE_{"vcfFile"};
  constexpr std::string_view VCF_FILE_GENOME_{"vcfGenome"};
  constexpr std::string_view VCF_INFO_EVIDENCE_{"vcfInfo"};
  constexpr std::string_view GENERAL_DATA_FILE_TYPE_{"generalFile"};

}


/// Parse data file definitions from the runtime XML configuration.
kgl::RuntimeDataFileMap kglx::parseDataFiles(const PropertyTree& tree, std::string_view work_directory) {

  RuntimeDataFileMap data_file_map;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not tree.getPropertyTreeVector(runtimeKey({DATA_FILE_LIST_}), property_tree_vector)) {

    ExecEnv::log().info("parseDataFiles; no Data files specified");
    return data_file_map;

  }

  for (auto const& [file_type_tag, sub_tree] : property_tree_vector) {

    if (file_type_tag == VCF_DATA_FILE_TYPE_) {

      XmlReader reader(sub_tree, "parseDataFiles(VCF)");

      std::string vcf_ident;
      std::string vcf_file_name;
      std::string vcf_parser_type;
      std::string vcf_reference_genome;
      std::string evidence_ident;
      reader.requiredString(DATA_FILE_IDENT_, vcf_ident)
            .requiredFile(DATA_FILE_NAME_, work_directory, vcf_file_name)
            .requiredString(DATA_PARSER_TYPE_, vcf_parser_type)
            .requiredString(VCF_FILE_GENOME_, vcf_reference_genome);

      if (not reader.ok()) {

        continue;

      }

      // Evidence is logged as missing but the VCF file object is still created.
      if (not reader.optionalString(VCF_INFO_EVIDENCE_, evidence_ident)) {

        ExecEnv::log().error("parseDataFiles; No VCF Info evidence specified for VCF file: {}", vcf_file_name);

      }

      auto file_info_ptr = std::make_shared<RuntimeVCFFileInfo>(vcf_ident, vcf_file_name,
                                                                vcf_parser_type, vcf_reference_genome,
                                                                evidence_ident);
      auto [iter, result] = data_file_map.try_emplace(vcf_ident, file_info_ptr);
      if (not result) {

        ExecEnv::log().error("parseDataFiles; Could not add VCF file ident: {} to map (duplicate)", vcf_ident);

      }

    } else if (file_type_tag == GENERAL_DATA_FILE_TYPE_) {

      XmlReader reader(sub_tree, "parseDataFiles(general)");

      std::string file_ident;
      std::string file_name;
      std::string parser_type;
      reader.requiredString(DATA_FILE_IDENT_, file_ident)
            .requiredFile(DATA_FILE_NAME_, work_directory, file_name)
            .requiredString(DATA_PARSER_TYPE_, parser_type);

      if (not reader.ok()) {

        continue;

      }

      auto file_info_ptr = std::make_shared<BaseFileInfo>(file_ident, file_name, parser_type);
      auto [iter, result] = data_file_map.try_emplace(file_ident, file_info_ptr);
      if (not result) {

        ExecEnv::log().error("parseDataFiles; Could not add file ident: {} to map (duplicate)", file_ident);

      }

    }

  }

  return data_file_map;

}