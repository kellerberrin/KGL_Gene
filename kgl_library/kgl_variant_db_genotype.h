//
// Created by kellerberrin on 1/04/18.
//

#ifndef KGL_VARIANT_DB_GENOTYPE_H
#define KGL_VARIANT_DB_GENOTYPE_H



#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant.h"
#include "kgl_variant_db_genome.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Genotype - A map of GenomeVariants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The GenomeVariant (haploid) map
using GenomeMap = std::map<GenomeId_t, std::shared_ptr<GenomeVariant>>;
class Genotype {

public:

  explicit Genotype(const GenotypeId_t& Genotype_id) : genotype_id_(Genotype_id) {}
  Genotype(const Genotype&) = delete; // Use deep copy.
  ~Genotype() = default;

  Genotype& operator=(const Genotype&) = delete; // Use deep copy.

  const GenotypeId_t& genotypeId() const { return genotype_id_; }
  void genotypeId(const GenotypeId_t& contig_id) { genotype_id_ = contig_id; }

  bool getGenome(const GenomeId_t& genome_id, std::shared_ptr<GenomeVariant>& genome) const;
  bool addGenome(std::shared_ptr<GenomeVariant> genome);

  const GenomeMap& getMap() const { return genome_map_; }


private:

  GenotypeId_t genotype_id_;
  GenomeMap genome_map_;

};


}   // namespace Genotype
}   // namespace kellerberrin



#endif //KGL_VARIANT_DB_GENOTYPE_H
