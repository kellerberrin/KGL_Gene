//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_GENOME_DB_H
#define KGL_GENOME_DB_H


#include "kgl_genome_genome.h"
#include "kgl_ontology_database.h"

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ResourceCollection - A map of different resources.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class Resource>
using ResourceMap = std::map<std::string, std::shared_ptr<const Resource>>;

template<class Resource>
class ResourceCollection {

public:

  explicit ResourceCollection() = default;
  ResourceCollection(const ResourceCollection&) = default;
  virtual ~ResourceCollection() = default;

  ResourceCollection& operator=(const ResourceCollection&) = default;

  // Returns false if the genome does not exist.
  [[nodiscard]] std::shared_ptr<const Resource> getResource(const std::string& resource_id) const {

    std::optional<std::shared_ptr<const Resource>> resource_opt = getOptionalResource(resource_id);
    if (not resource_opt) {

      ExecEnv::log().critical("ResourceCollection::getOptionalResource; resource: {} not found", resource_id);

    }

    return resource_opt.value();

  }

  [[nodiscard]] std::optional<std::shared_ptr<const Resource>> getOptionalResource(const std::string& resource_id) const {

    auto result = resource_map_.find(resource_id);
    if (result != resource_map_.end()) {

      return result->second;

    } else {

      return std::nullopt;

    }

  }

  [[nodiscard]] const ResourceMap<Resource>& getMap() const { return resource_map_; }

  // Returns false if the genome already exists.
  [[nodiscard]] bool addResource(const std::string& resource_id, const std::shared_ptr<const Resource>& resource_ptr) {

    if (not resource_ptr or resource_id.empty()) {

      ExecEnv::log().error("ResourceCollection::addResource; attempt to add resource: {} with nullptr or empty resource id.", resource_id);
      return false;

    }

    auto [it, result] = resource_map_.try_emplace(resource_id, resource_ptr);

    return result;

  }

private:

  // A map of all active genome databases.
  ResourceMap<Resource> resource_map_;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeCollection - A map of different organism genomes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeCollection = ResourceCollection<GenomeReference>;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OntologyCollection - A map of different Ontology databases.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using OntologyDatabase = kellerberrin::ontology::OntologyDatabase;
using OntologyCollection = ResourceCollection<OntologyDatabase>;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Package Resources. A container object to hold the resources specified by an analysis package.
// This object is presented to the analysis package on initialization.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AnalysisResources {

public:

  AnalysisResources() = default;
  ~AnalysisResources() = default;

// Genome resources

  [[nodiscard]] const GenomeCollection& getGenomes() const { return genome_resources_; }
  [[nodiscard]] bool addGenome(const std::shared_ptr<const GenomeReference>& genome_database) {

    return genome_resources_.addResource(genome_database->genomeId(), genome_database);

  }

// Gene Ontology Resources

  [[nodiscard]] const OntologyCollection& getOntologies() const { return ontology_resources_; }
  [[nodiscard]] bool addOntology(std::shared_ptr<const OntologyDatabase>& ontology_database) {

    return ontology_resources_.addResource(ontology_database->ontologyIdent(), ontology_database);

  }


private:

  GenomeCollection genome_resources_;
  OntologyCollection ontology_resources_;

};



}   // end namespace

#endif //KGL_GENOME_DB_H
