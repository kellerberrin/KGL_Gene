//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_xml_parsers.h"
#include "kgl_xml_reader.h"
#include "kgl_xml_detail.h"
#include "kel_utility.h"

#include <set>
#include <utility>


namespace kgl = kellerberrin::genome;
namespace kglx = kellerberrin::genome::xml;

namespace {

  // Active package tags.
  constexpr std::string_view EXECUTE_LIST_{"executeList"};
  constexpr std::string_view PACKAGE_{"package"};

  // Contig alias tags.
  constexpr std::string_view ALIAS_LIST_{"aliasList"};
  constexpr std::string_view ALIAS_IDENT_{"ident"};
  constexpr std::string_view ALIAS_TYPE_{"chromosomeType"};
  constexpr std::string_view ALIAS_ENTRY_{"alias"};

  // VCF evidence tags.
  constexpr std::string_view EVIDENCE_LIST_{"evidenceList"};
  constexpr std::string_view EVIDENCE_IDENT_{"evidenceIdent"};
  constexpr std::string_view EVIDENCE_INFO_LIST_{"vcfInfoList"};
  constexpr std::string_view EVIDENCE_INFO_ITEM_{"vcfInfoItem"};

}


/// Parse the active package list from the runtime XML configuration.
kgl::ActivePackageVector kglx::parseActivePackages(const PropertyTree& tree) {

  ActivePackageVector active_packages;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not tree.getPropertyTreeVector(runtimeKey({EXECUTE_LIST_, ACTIVE_}), property_tree_vector)) {

    ExecEnv::log().info("parseActivePackages; No Active Packages Specified");
    return active_packages;

  }

  for (auto const& [tag, sub_tree] : property_tree_vector) {

    if (tag == PACKAGE_) {

      active_packages.emplace_back(sub_tree.getValue());

    }

  }

  return active_packages;

}


/// Parse contig alias mappings from the runtime XML configuration.
kgl::ContigAliasMap kglx::parseContigAliases(const PropertyTree& tree) {

  ContigAliasMap contig_alias_map;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not tree.getPropertyTreeVector(runtimeKey({ALIAS_LIST_}), property_tree_vector)) {

    ExecEnv::log().info("parseContigAliases; No Contig Alias Specified");
    return contig_alias_map;

  }

  for (auto const& [sub_tree_id, sub_tree] : property_tree_vector) {

    XmlReader reader(sub_tree, "parseContigAliases");

    std::string contig_ident;
    std::string chromosome_type;
    reader.requiredString(ALIAS_IDENT_, contig_ident)
          .requiredString(ALIAS_TYPE_, chromosome_type);
    if (not reader.ok()) {

      continue;

    }

    contig_alias_map.setAlias(contig_ident, contig_ident, chromosome_type);

    std::vector<std::string> alias_vector;
    if (sub_tree.checkProperty(std::string(ALIAS_ENTRY_))) {

      if (not sub_tree.getNodeVector(std::string(ALIAS_ENTRY_), alias_vector)) {

        ExecEnv::log().warn("parseContigAliases; No Alias Specified for Contig: {}", sub_tree_id);

      }

    }

    for (auto const& alias : alias_vector) {

      contig_alias_map.setAlias(alias, contig_ident, chromosome_type);

    }

  }

  return contig_alias_map;

}


/// Parse VCF evidence mappings from the runtime XML configuration.
kgl::VariantEvidenceMap kglx::parseEvidenceMap(const PropertyTree& tree) {

  VariantEvidenceMap variant_evidence_map;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not tree.getPropertyTreeVector(runtimeKey({EVIDENCE_LIST_}), property_tree_vector)) {

    ExecEnv::log().info("parseEvidenceMap; No Variant Evidence items specified.");
    return variant_evidence_map;

  }

  for (auto const& [tag, sub_tree] : property_tree_vector) {

    XmlReader reader(sub_tree, "parseEvidenceMap");

    std::string evidence_ident;
    reader.requiredString(EVIDENCE_IDENT_, evidence_ident);
    if (not reader.ok()) {

      continue;

    }

    std::vector<SubPropertyTree> info_tree_vector;
    if (not sub_tree.getPropertyTreeVector(std::string(EVIDENCE_INFO_LIST_), info_tree_vector)) {

      ExecEnv::log().info("parseEvidenceMap; no info items specified for evidence ident: {}", evidence_ident);

    }

    std::set<std::string> evidence_list;
    for (auto const& [info_tag, info_sub_tree] : info_tree_vector) {

      if (info_tag != EVIDENCE_INFO_ITEM_) continue;

      auto [iter, inserted] = evidence_list.insert(info_sub_tree.getValue());
      if (not inserted) {

        ExecEnv::log().warn("parseEvidenceMap; Duplicate Info item: {} for evidence: {}",
                            info_sub_tree.getValue(), evidence_ident);

      }

    }

    variant_evidence_map.setEvidence(evidence_ident, evidence_list);

  }

  return variant_evidence_map;

}