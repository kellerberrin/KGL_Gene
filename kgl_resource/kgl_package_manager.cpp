//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_package_manager.h"
#include "kgl_resource_factory.h"
#include "kgl_variant_factory_parsers.h"


namespace kgl = kellerberrin::genome;


/// Load all resources requested by a package using the ResourceFactory registry.
std::shared_ptr<const kgl::AnalysisResources> kgl::PackageManager::loadResources(const RuntimePackage& package) const {

  auto resource_ptr = std::make_shared<AnalysisResources>();

  for (auto const& [resource_type, resource_ident] : package.resourceList()) {

    auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, resource_ident);
    if (not params_opt) {

      ExecEnv::log().critical("PackageManager::loadResources; Package: {}, Resource Type: '{}', Ident: '{}' not defined",
                               package.packageIdentifier(), resource_type, resource_ident);

    }

    const ResourceFactory* factory = findResourceFactory(resource_type);
    if (not factory) {

      ExecEnv::log().critical("PackageManager::loadResources; Package: {}, no factory for resource type: '{}'",
                              package.packageIdentifier(), resource_type);

    }

    resource_ptr->addResource(factory->load(params_opt.value()));

  }

  return resource_ptr;

}


/// Read a single data file identified by its package iteration list.
std::shared_ptr<kgl::DataDB>
kgl::PackageManager::readDataFile(const RuntimePackage& package,
                                  const std::shared_ptr<const AnalysisResources>& resource_ptr,
                                  const std::string& data_file_ident) const {

  ExecEnv::log().info("Package: {}, Data file ident: {}", package.packageIdentifier(), data_file_ident);

  auto result = runtime_config_.dataFileMap().find(data_file_ident);
  if (result == runtime_config_.dataFileMap().end()) {

    ExecEnv::log().critical("PackageManager::readDataFile; Package: {}, data file ident: {}, not defined",
                             package.packageIdentifier(), data_file_ident);

  }

  auto const& [file_ident, file_info_ptr] = *result;

  return ParserSelection::parseData(resource_ptr, file_info_ptr,
                                     runtime_config_.evidenceMap(), runtime_config_.contigAlias());

}