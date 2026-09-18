//
// Created by kellerberrin on 26/08/26.
//

#ifndef KGL_kgl_resource_pf3kcoi_H
#define KGL_kgl_resource_pf3kcoi_H


#include "kgl_resource_factory.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/// Factory for Pf3K COI (complexity of infection) data resources (Pf3KCOIFile).
class Pf3KCOIResourceFactory : public ResourceFactory {

public:

  [[nodiscard]] std::string_view resourceType() const override;
  [[nodiscard]] std::optional<ResourceParameters> parseXML(const PropertyTree& sub_tree,
                                                          std::string_view work_directory) const override;
  [[nodiscard]] std::shared_ptr<const ResourceBase> load(const ResourceParameters& parameters) const override;

};


}   // end namespace


#endif //KGL_kgl_resource_pf3kcoi_H
