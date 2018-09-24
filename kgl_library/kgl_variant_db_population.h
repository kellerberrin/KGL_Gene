//
// Created by kellerberrin on 8/01/18.
//

#ifndef KGL_VARIANT_DB_POPULATION_H
#define KGL_VARIANT_DB_POPULATION_H




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
// Simple container to hold phased genome variants for populations
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using PopulationVariantMap = std::map<GenomeId_t, std::shared_ptr<const GenomeVariant>>;
class PhasedPopulation {

public:

  explicit PhasedPopulation(const PopulationId_t& population_id) : population_id_(population_id) {}
  PhasedPopulation(const PhasedPopulation&) = delete;
  virtual ~PhasedPopulation() = default;

  PhasedPopulation& operator=(const PhasedPopulation&) = delete; // Use deep copy.

  // Create the genome variant if it does not exist.
  bool getCreateGenome(const GenomeId_t& genome_id,
                       PhaseId_t ploidy,
                       const std::shared_ptr<const GenomeDatabase>& genome_db,
                       std::shared_ptr<GenomeVariant>& genome);

  // Returns false if the genome does not exist.
  bool getGenomeVariant(const GenomeId_t& genome_id, std::shared_ptr<const GenomeVariant>& genome_variant) const;

  bool addGenomeVariant(std::shared_ptr<const GenomeVariant> genome_variant);

  size_t variantCount() const;

  const PopulationVariantMap& getMap() const { return population_variant_map_; }

  std::shared_ptr<PhasedPopulation> filterVariants(const VariantFilter& filter) const;

  std::shared_ptr<PhasedPopulation> filterGenomes(const PopulationId_t& population_id, const std::vector<GenomeId_t>& list) const;

  const PopulationId_t& populationId() const { return population_id_; }

private:

  PopulationVariantMap population_variant_map_;
  PopulationId_t population_id_;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_DB_POPULATION_H
