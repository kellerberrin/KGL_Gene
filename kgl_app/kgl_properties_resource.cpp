//
// Created by kellerberrin on 31/03/23.
//


#include "kgl_properties_resource.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;


// A map of resources.
kgl::ResourceDefinitions kgl::ResourceProperties::getRuntimeResources() const {

  ResourceDefinitions resource_definitions;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + std::string(RESOURCE_LIST_);
  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_ptr_->getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().warn("RuntimeProperties::getRuntimeResources; No runtime resources specified.");
    return resource_definitions; // return empty map.

  }

  for (const auto &[tree_type, sub_tree]: property_tree_vector) {

    std::optional<ResourceParameters> resource_opt{std::nullopt};
    switch (Utility::hash(tree_type)) {

      case Utility::hash(GENOME_RESOURCE_ID_):
        resource_opt = genomeResourceXML(sub_tree);
        break;

      case Utility::hash(ONTOLOGY_RESOURCE_ID_):
        resource_opt = ontologyResourceXML(sub_tree);
        break;

      case Utility::hash(GENE_NOMENCLATURE_RESOURCE_ID_):
        resource_opt = geneIDResourceXML(sub_tree);
        break;

      case Utility::hash(GENEALOGY_RESOURCE_ID_):
        resource_opt = genealogyIDResourceXML(sub_tree);
        break;

      case Utility::hash(CITATION_RESOURCE_ID_):
        resource_opt = citationResourceXML(sub_tree);
        break;

      case Utility::hash(ENTREZ_RESOURCE_ID_):
        resource_opt = entrezResourceXML(sub_tree);
        break;

      case Utility::hash(PF7SAMPLE_RESOURCE_ID_):
        resource_opt = Pf7SampleResourceXML(sub_tree);
        break;

      case Utility::hash(PF7FWS_RESOURCE_ID_):
        resource_opt = Pf7FwsResourceXML(sub_tree);
        break;

      case Utility::hash(PF7DISTANCE_RESOURCE_ID_):
        resource_opt = Pf7DistanceResourceXML(sub_tree);
        break;

      case Utility::hash(PUBMED_API_RESOURCE_ID_):
        resource_opt = PubmedLitAPIResourceXML(sub_tree);
        break;

      case Utility::hash(GENOMEAUX_RESOURCE_ID_):
        resource_opt = auxIDResourceXML(sub_tree);
        break;

      case Utility::hash(PF3K_COI_RESOURCE_ID_):
        resource_opt = Pf3KCOIResourceXML(sub_tree);
        break;

      default:
        ExecEnv::log().warn("RuntimeProperties::getRuntimeResources, Unexpected resource: '{}'; ignored", tree_type);
        break;

    }

    if (resource_opt) {

      resource_definitions.insert({tree_type, resource_opt.value()});

    }

  }

  return resource_definitions;

}

