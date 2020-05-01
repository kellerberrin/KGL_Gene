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


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container to hold phased genome variants for populations
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using PopulationVariantMap = std::map<GenomeId_t, std::shared_ptr<const GenomeVariant>>;
class PhasedPopulation {

public:

  explicit PhasedPopulation(const PopulationId_t& population_id) : population_id_(population_id) {}
  PhasedPopulation(const PhasedPopulation&) = delete;
  virtual ~PhasedPopulation() = default;

  [[nodiscard]] PhasedPopulation& operator=(const PhasedPopulation&) = delete; // Use deep copy.

  // Create the genome variant if it does not exist.
  [[nodiscard]] bool getCreateGenome( const GenomeId_t& genome_id,
                                      PhaseId_t ploidy,
                                      const std::shared_ptr<const GenomeReference>& genome_db,
                                      std::shared_ptr<GenomeVariant>& genome);

  // Returns false if the genome does not exist.
  [[nodiscard]] bool getGenomeVariant(const GenomeId_t& genome_id, std::shared_ptr<const GenomeVariant>& genome_variant) const;

  [[nodiscard]] bool addGenomeVariant(std::shared_ptr<const GenomeVariant> genome_variant);

  [[nodiscard]] size_t variantCount() const;

  [[nodiscard]] const PopulationVariantMap& getMap() const { return population_variant_map_; }

  [[nodiscard]] std::shared_ptr<PhasedPopulation> filterVariants(const VariantFilter& filter) const;

  [[nodiscard]] std::shared_ptr<PhasedPopulation> filterGenomes(const PopulationId_t& population_id, const std::vector<GenomeId_t>& list) const;

  // First element of the pair is the population genome name which is renamed to the second element of the pair
  [[nodiscard]] std::shared_ptr<PhasedPopulation> filterRenameGenomes( const PopulationId_t& population_id,
                                                                       const std::vector<std::pair<GenomeId_t, GenomeId_t >>& source_pairs) const;


  [[nodiscard]] const PopulationId_t& populationId() const { return population_id_; }

private:

  PopulationVariantMap population_variant_map_;
  PopulationId_t population_id_;

};



}   // end namespace



#endif //KGL_VARIANT_DB_POPULATION_H
