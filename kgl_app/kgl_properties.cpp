//
// Created by kellerberrin on 11/11/18.
//

#include "kgl_properties.h"
#include "kgl_xml_parsers.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;


/// Reads the runtime properties from the specified XML file.
bool kgl::RuntimeProperties::readProperties(const std::string& properties_file) {

  std::string properties_path = Utility::filePath(properties_file, work_directory_);
  return property_tree_ptr_->readProperties(properties_path);

}


/// Returns the vector of active packages parsed from the XML properties.
const kgl::ActivePackageVector& kgl::RuntimeProperties::getActivePackages() const {

  return cache(active_packages_cache_, [&] { return xml::parseActivePackages(*property_tree_ptr_); });

}


/// Returns the map of runtime packages parsed from the XML properties.
const kgl::RuntimePackageMap& kgl::RuntimeProperties::getPackageMap() const {

  return cache(package_map_cache_, [&] { return xml::parsePackageMap(*property_tree_ptr_); });

}


/// Returns the map of runtime analysis definitions parsed from the XML properties.
const kgl::RuntimeAnalysisMap& kgl::RuntimeProperties::getAnalysisMap() const {

  return cache(analysis_map_cache_, [&] { return xml::parseAnalysisMap(*property_tree_ptr_); });

}


/// Returns the resource definitions parsed from the XML properties.
const kgl::ResourceDefinitions& kgl::RuntimeProperties::getRuntimeResources() const {

  return cache(runtime_resources_cache_, [&] { return xml::parseResources(*property_tree_ptr_, work_directory_); });

}


/// Returns the map of data files parsed from the XML properties.
const kgl::RuntimeDataFileMap& kgl::RuntimeProperties::getDataFiles() const {

  return cache(data_files_cache_, [&] { return xml::parseDataFiles(*property_tree_ptr_, work_directory_); });

}


/// Returns the contig alias map parsed from the XML properties.
const kgl::ContigAliasMap& kgl::RuntimeProperties::getContigAlias() const {

  return cache(contig_alias_cache_, [&] { return xml::parseContigAliases(*property_tree_ptr_); });

}


/// Returns the variant evidence map parsed from the XML properties.
const kgl::VariantEvidenceMap& kgl::RuntimeProperties::getEvidenceMap() const {

  return cache(evidence_map_cache_, [&] { return xml::parseEvidenceMap(*property_tree_ptr_); });

}


/// Returns the active parameter list parsed from the XML properties.
const kgl::ActiveParameterList& kgl::RuntimeProperties::getParameterMap() const {

  return cache(parameter_map_cache_, [&] { return xml::parseParameters(*property_tree_ptr_); });

}