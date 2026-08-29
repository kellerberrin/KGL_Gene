//
// Created by kellerberrin on 26/08/26.
//

#ifndef KGL_kgl_resource_aux_H
#define KGL_kgl_resource_aux_H


#include "kgl_resource_factory.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/// Factory for genome auxiliary information resources (auxFile).
class AuxResourceFactory : public ResourceFactory {

public:

  [[nodiscard]] std::string_view resourceType() const override;
  [[nodiscard]] std::optional<ResourceParameters> parseXML(const PropertyTree& sub_tree,
                                                          std::string_view work_directory) const override;
  [[nodiscard]] std::shared_ptr<const ResourceBase> load(const ResourceParameters& parameters) const override;

};


}   // end namespace


#endif //KGL_kgl_resource_aux_H
