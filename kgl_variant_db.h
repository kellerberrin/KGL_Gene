//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_VARIANT_DB_H
#define KGL_VARIANT_DB_H


#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant.h"
#include "kgl_variant_db_contig.h"
#include "kgl_variant_db_genome.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container to hold genome variants for populations
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using PopulationVariantMap = std::map<GenomeId_t, std::shared_ptr<const GenomeVariant>>;
class PopulationVariant {

public:

  explicit PopulationVariant(const std::string& population_id) : population_id_(population_id) {}
  PopulationVariant(const PopulationVariant&) = default;
  ~PopulationVariant() = default;

  bool addGenomeVariant(std::shared_ptr<const GenomeVariant> genome_variant);

  bool getGenomeVariant(const GenomeId_t& genome_id, std::shared_ptr<const GenomeVariant>& genome_variant);

  const PopulationVariantMap& getMap() const { return population_variant_map_; }

private:

  PopulationVariantMap population_variant_map_;
  std::string population_id_;

};


}   // namespace genome
}   // namespace kellerberrin


// Not in kgl:: namespace.
std::ostream & operator<<(std::ostream &os, const kellerberrin::genome::GenomeVariant& genome_variant);


#endif //KGL_VARIANT_DB_H
