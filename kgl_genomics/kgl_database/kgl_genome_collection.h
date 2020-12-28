//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_GENOME_DB_H
#define KGL_GENOME_DB_H


#include "kgl_genome_genome.h"

#include <memory>
#include <string>
#include <vector>
#include <map>


namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeCollection - A map of different organism genomes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeMap = std::map<GenomeId_t, std::shared_ptr<const GenomeReference>>;

class GenomeCollection {

public:

  explicit GenomeCollection() = default;
  GenomeCollection(const GenomeCollection&) = default;
  virtual ~GenomeCollection() = default;

  GenomeCollection& operator=(const GenomeCollection&) = default;

  // High level function creates a collection of genomes.
  [[nodiscard]] static std::shared_ptr<GenomeCollection> createGenomeCollection(const RuntimeProperties& runtime_options);

  // Returns false if the genome does not exist.
  [[nodiscard]] std::shared_ptr<const GenomeReference> getGenome(const std::string& GenomeID) const;
  [[nodiscard]] std::optional<std::shared_ptr<const GenomeReference>> getOptionalGenome(const GenomeId_t& genome_id) const;

  [[nodiscard]] const GenomeMap& getMap() const { return genome_map_; }

  // Returns false if the genome already exists.
  [[nodiscard]] bool addGenome(std::shared_ptr<const GenomeReference> genome_database);

private:

  // A map of all active genome databases.
  GenomeMap genome_map_;


};


}   // end namespace


#endif //KGL_GENOME_DB_H
