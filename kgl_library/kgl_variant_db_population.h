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
#include "kgl_variant_db_genotype.h"


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
  virtual ~PopulationVariant() = default;

  // Returns false if the genome does not exist.
  bool getGenomeVariant(const GenomeId_t& genome_id, std::shared_ptr<const GenomeVariant>& genome_variant) const;

  bool addGenomeVariant(std::shared_ptr<const GenomeVariant> genome_variant);

  const PopulationVariantMap& getMap() const { return population_variant_map_; }

private:

  PopulationVariantMap population_variant_map_;
  std::string population_id_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container to hold genomes with a ploidy > 1
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using GenotypeMap = std::map<GenomeId_t, std::shared_ptr<Genotype>>;
class GenotypePopulation {

public:

  explicit GenotypePopulation(const std::string& population_id) : population_id_(population_id) {}
  GenotypePopulation(const GenotypePopulation&) = default;
  virtual ~GenotypePopulation() = default;

  // Returns false if the genome does not exist.
  bool getGenotype(const GenotypeId_t& genotype_id, std::shared_ptr<Genotype>& genotype) const;
  bool addGenotype(std::shared_ptr<Genotype> genotype);

  const GenotypeMap& getMap() const { return genotype_map_; }

private:

  GenotypeMap genotype_map_;
  std::string population_id_;

};


}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_DB_POPULATION_H
