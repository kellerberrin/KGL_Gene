//
// Created by kellerberrin on 13/08/18.
//

#ifndef KGL_VARIANT_DB_UNPHASED_POPULATION_H
#define KGL_VARIANT_DB_UNPHASED_POPULATION_H


#include "kgl_variant_db_unphased.h"


namespace kellerberrin::genome {   //  organization::project


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants until they can be phased.
// This object hold variants for a population.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


using UnphasedGenomeMap = std::map<GenomeId_t, std::shared_ptr<UnphasedGenome>>;
class UnphasedPopulation {

public:

  explicit UnphasedPopulation() = default;
  UnphasedPopulation(const UnphasedPopulation&) = delete; // Use deep copy.
  virtual ~UnphasedPopulation() = default;

  UnphasedPopulation& operator=(const UnphasedPopulation&) = delete; // Use deep copy.

  // Use this to copy the object.
  [[nodiscard]] std::shared_ptr<UnphasedPopulation> deepCopy() const;

  // Create the genome variant if it does not exist.
  [[nodiscard]] bool getCreateGenome(const GenomeId_t& genome_id, std::shared_ptr<UnphasedGenome>& genome);

  [[nodiscard]] size_t variantCount() const;
  void popStatistics() const; // output to logger
  [[nodiscard]] std::vector<GenomeId_t> genomeList() const;

  [[nodiscard]] std::shared_ptr<UnphasedPopulation> filterVariants(const VariantFilter& filter) const;

  [[nodiscard]] const UnphasedGenomeMap& getMap() const { return genome_map_; }

  void clear() { genome_map_.clear(); }

  // Generate phasing statitics.
  [[nodiscard]] bool genomePhasingStats( const GenomeId_t& genome_id,
                                         size_t& heterozygous,
                                         size_t& homozygous,
                                         size_t& singleheterozygous) const;

  [[nodiscard]] bool addGenome(std::shared_ptr<UnphasedGenome> genome);

private:

  UnphasedGenomeMap genome_map_;


};


}   // end namespace

[[nodiscard]] std::ostream& operator<<(std::ostream& ostream, std::shared_ptr<const kellerberrin::genome::UnphasedPopulation> unphased_ptr);
[[nodiscard]] std::ostream& operator<<(std::ostream& ostream, const kellerberrin::genome::UnphasedPopulation& unphased);


#endif //KGL_VARIANT_DB_UNPHASED_POPULATION_H
