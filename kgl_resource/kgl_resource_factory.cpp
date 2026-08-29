//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_factory.h"
#include "kgl_resource_type.h"
#include "kgl_resource_genome.h"
#include "kgl_resource_ontology.h"
#include "kgl_resource_gene_nomenclature.h"
#include "kgl_resource_genealogy.h"
#include "kgl_resource_aux.h"
#include "kgl_resource_citation.h"
#include "kgl_resource_entrez.h"
#include "kgl_resource_pubmed.h"
#include "kgl_resource_pf7_sample.h"
#include "kgl_resource_pf7_fws.h"
#include "kgl_resource_pf7_distance.h"
#include "kgl_resource_pf3kcoi.h"

#include "kel_exec_env.h"

namespace kgl = kellerberrin::genome;


/// Locate the factory for a given resource type. Returns nullptr for unknown types.
const kgl::ResourceFactory* kgl::findResourceFactory(std::string_view resource_type) {

  // Function-local static; initialized once on first access (thread-safe since C++11).
  // std::less<> enables transparent string_view lookup without allocating a std::string key.
  static const std::map<std::string_view, std::unique_ptr<const ResourceFactory>, std::less<>> factories = [] {

    std::map<std::string_view, std::unique_ptr<const ResourceFactory>, std::less<>> m;
    m.emplace(ResourceType::GENOME,            std::make_unique<GenomeResourceFactory>());
    m.emplace(ResourceType::ONTOLOGY,          std::make_unique<OntologyResourceFactory>());
    m.emplace(ResourceType::GENE_NOMENCLATURE, std::make_unique<GeneNomenclatureFactory>());
    m.emplace(ResourceType::GENEALOGY,        std::make_unique<GenealogyResourceFactory>());
    m.emplace(ResourceType::GENOME_AUX,       std::make_unique<AuxResourceFactory>());
    m.emplace(ResourceType::CITATION,         std::make_unique<CitationResourceFactory>());
    m.emplace(ResourceType::ENTREZ,           std::make_unique<EntrezResourceFactory>());
    m.emplace(ResourceType::PUBMED_API,       std::make_unique<PubmedResourceFactory>());
    m.emplace(ResourceType::PF7_SAMPLE,       std::make_unique<Pf7SampleResourceFactory>());
    m.emplace(ResourceType::PF7_FWS,          std::make_unique<Pf7FwsResourceFactory>());
    m.emplace(ResourceType::PF7_DISTANCE,     std::make_unique<Pf7DistanceResourceFactory>());
    m.emplace(ResourceType::PF3K_COI,         std::make_unique<Pf3KCOIResourceFactory>());
    return m;

  }();

  auto it = factories.find(resource_type);
  if (it == factories.end()) {

    return nullptr;

  }

  return it->second.get();

}