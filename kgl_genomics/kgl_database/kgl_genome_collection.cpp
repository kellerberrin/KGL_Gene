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

  auto result = genome_map_.insert(std::pair<GenomeId_t, std::shared_ptr<const GenomeReference>>(genome_database->genomeId(), genome_database));

  return result.second;

}


std::shared_ptr<const kgl::GenomeReference> kgl::GenomeCollection::getGenome(const std::string& GenomeID) const {

  std::optional<std::shared_ptr<const GenomeReference>> genome_opt = getOptionalGenome(GenomeID);
  if (not genome_opt) {

    ExecEnv::log().critical("GenomeCollection::getOptionalGenome; genome: {} not found", GenomeID);

  }

  return genome_opt.value();

}



std::shared_ptr<kgl::GenomeCollection> kgl::GenomeCollection::createGenomeCollection(const RuntimeProperties& runtime_options) {

  std::vector<std::string> genome_list;
  if (not runtime_options.getActiveGenomes(genome_list)) {

    ExecEnv::log().error("GenomeCollection::createGenomeCollection; Problem retrieving runtime genome list");

  }

  std::shared_ptr<kgl::GenomeCollection> genome_collection(std::make_shared<kgl::GenomeCollection>());

  for (auto genome : genome_list) {

    // Create the genome database.
    std::shared_ptr<GenomeReference> genome_ptr = GenomeReference::createGenomeDatabase(runtime_options, genome);

    if (not genome_collection->addGenome(genome_ptr)) {

      ExecEnv::log().error("GenomeCollection::createGenomeCollection; Unable to add Genome Database: {} (probable duplicate)", genome_ptr->genomeId());

    }

  }

  return genome_collection;

}

