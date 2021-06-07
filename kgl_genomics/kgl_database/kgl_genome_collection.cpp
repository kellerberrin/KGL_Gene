//
// Created by kellerberrin on 7/6/21.
//

#include "kgl_genome_collection.h"

namespace kgl = kellerberrin::genome;

// Returns false if the genome does not exist.
[[nodiscard]] std::shared_ptr<const kgl::GenomeReference> kgl::GenomeCollection::getGenome(const std::string &resource_id) const {

  std::optional<std::shared_ptr<const GenomeReference>> resource_opt = getOptionalGenome(resource_id);
  if (not resource_opt) {

    ExecEnv::log().critical("GenomeCollection::getOptionalGenome; resource: {} not found", resource_id);

  }

  return resource_opt.value();

}



[[nodiscard]] std::optional<std::shared_ptr<const kgl::GenomeReference>> kgl::GenomeCollection::getOptionalGenome(const std::string &resource_id) const {

  auto result = genome_map_.find(resource_id);
  if (result != genome_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}


// Returns false if the genome already exists.
bool kgl::GenomeCollection::addGenome(const std::shared_ptr<const GenomeReference> &genome_ptr) {

  if (not genome_ptr) {

    ExecEnv::log().error("ResourceCollection::addGenome; attempt to add resource: {} with nullptr or empty resource id.", genome_ptr->genomeId());
    return false;

  }

  auto[it, result] = genome_map_.try_emplace(genome_ptr->genomeId(), genome_ptr);

  return result;

}
