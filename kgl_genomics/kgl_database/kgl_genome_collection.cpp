//
// Created by kellerberrin on 28/12/20.
//

#include "kgl_genome_collection.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeCollection - A map of different organism genomes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::optional<std::shared_ptr<const kgl::GenomeReference>> kgl::GenomeCollection::getOptionalGenome(const GenomeId_t& genome_id) const {

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}


bool kgl::GenomeCollection::addGenome(std::shared_ptr<const GenomeReference> genome_database) {

  if (not genome_database) {

    ExecEnv::log().error("GenomeCollection::addGenome; attempt to add genome with nullptr");
    return false;

  }

  auto [it, result] = genome_map_.try_emplace(genome_database->genomeId(), genome_database);

  return result;

}


std::shared_ptr<const kgl::GenomeReference> kgl::GenomeCollection::getGenome(const std::string& GenomeID) const {

  std::optional<std::shared_ptr<const GenomeReference>> genome_opt = getOptionalGenome(GenomeID);
  if (not genome_opt) {

    ExecEnv::log().critical("GenomeCollection::getOptionalGenome; genome: {} not found", GenomeID);

  }

  return genome_opt.value();

}


