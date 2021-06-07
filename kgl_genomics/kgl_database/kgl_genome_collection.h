//
// Created by kellerberrin on 7/6/21.
//

#ifndef KGL_GENOME_COLLECTION_H
#define KGL_GENOME_COLLECTION_H

#include "kgl_genome_genome.h"
#include "kgl_ontology_database.h"

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeCollection - A map of different genomes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using GenomeReferenceMap = std::map<std::string, std::shared_ptr<const GenomeReference>>;

class GenomeCollection {

public:

  GenomeCollection() = default;
  GenomeCollection(const GenomeCollection &) = default;
  virtual ~GenomeCollection() = default;

  GenomeCollection &operator=(const GenomeCollection &) = default;

  [[nodiscard]] std::shared_ptr<const GenomeReference> getGenome(const std::string &genome_id) const;
  [[nodiscard]] std::optional<std::shared_ptr<const GenomeReference>> getOptionalGenome(const std::string &genome_id) const;
  [[nodiscard]] const GenomeReferenceMap &getMap() const { return genome_map_; }
  // Returns false if the genome already exists.
  bool addGenome(const std::shared_ptr<const GenomeReference> &genome_ptr);

private:

  // A map of all active genome databases.
  GenomeReferenceMap genome_map_;


};


}


#endif //KGL_GENOME_COLLECTION_H
