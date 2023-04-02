//
// Created by kellerberrin on 31/03/23.
//


#include "kgl_properties_resource.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;


// A map of resources.
kgl::ResourceDefinitions kgl::ResourceProperties::getRuntimeResources() const {

  ResourceDefinitions resource_definitions;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + RESOURCE_LIST_;
  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_ptr_->getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().warn("RuntimeProperties::getRuntimeResources; No runtime resources specified.");
    return resource_definitions; // return empty map.

  }

  for (const auto &[tree_type, sub_tree]: property_tree_vector) {

    std::optional<ResourceParameters> resource_opt{std::nullopt};
    switch (Utility::hash(tree_type)) {

      case Utility::hash(GENOME_DATABASE_):
        resource_opt = genomeResource(sub_tree);
        break;

      case Utility::hash(ONTOLOGY_DATABASE_):
        resource_opt = ontologyDatabase(sub_tree);
        break;

      case Utility::hash(GENE_ID_DATABASE_):
        resource_opt = geneIDDatabase(sub_tree);
        break;

      case Utility::hash(GENEALOGY_ID_DATABASE_):
        resource_opt = genealogyIDDatabase(sub_tree);
        break;

      case Utility::hash(CITATION_DATABASE_):
        resource_opt = citationDatabase(sub_tree);
        break;

      case Utility::hash(ENTREZ_DATABASE_):
        resource_opt = entrezDatabase(sub_tree);
        break;

      case Utility::hash(PF7_SAMPLE_DATABASE_):
        resource_opt = Pf7SampleDatabase(sub_tree);
        break;

      case Utility::hash(PUBMED_LIT_API_):
        resource_opt = PubmedLitAPI(sub_tree);
        break;

      case Utility::hash(AUX_ID_DATABASE_):
        resource_opt = auxIDDatabase(sub_tree);
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

