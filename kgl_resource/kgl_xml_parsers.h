//
// Created by kellerberrin on 26/08/26.
//

#ifndef KGL_XML_PARSERS_H
#define KGL_XML_PARSERS_H


#include "kgl_runtime.h"
#include "kgl_runtime_resource.h"
#include "kel_property_tree.h"

#include <string_view>


namespace kellerberrin::genome::xml {   //  organization::project::xml level namespace


/// XML section parser free functions. Each parses one section of the runtime XML
/// configuration and returns the corresponding data structure.

[[nodiscard]] ActivePackageVector parseActivePackages(const PropertyTree& tree);
[[nodiscard]] RuntimePackageMap parsePackageMap(const PropertyTree& tree);
[[nodiscard]] RuntimeAnalysisMap parseAnalysisMap(const PropertyTree& tree);
[[nodiscard]] RuntimeDataFileMap parseDataFiles(const PropertyTree& tree, std::string_view work_directory);
[[nodiscard]] ContigAliasMap parseContigAliases(const PropertyTree& tree);
[[nodiscard]] VariantEvidenceMap parseEvidenceMap(const PropertyTree& tree);
[[nodiscard]] ActiveParameterList parseParameters(const PropertyTree& tree);
[[nodiscard]] ResourceDefinitions parseResources(const PropertyTree& tree, std::string_view work_directory);


}   // end namespace


#endif //KGL_XML_PARSERS_H