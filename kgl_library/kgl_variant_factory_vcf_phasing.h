//
// Created by kellerberrin on 2/04/18.
//

#ifndef KGL_VARIANT_FACTORY_VCF_PHASING_H
#define KGL_VARIANT_FACTORY_VCF_PHASING_H


#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds multi ploid variants until they can be phased.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using ContigMultiMap = std::multimap<ContigOffset_t, std::shared_ptr<Variant>>;
using VCFContigMap = std::map<ContigId_t, std::shared_ptr<ContigMultiMap>>;
class VCFGenome {

public:

  explicit VCFGenome() = default;
  VCFGenome(const VCFGenome&) = default;
  virtual ~VCFGenome() = default;

  size_t variantCount() const;

  bool addVariant(std::shared_ptr<Variant> variant);

  const VCFContigMap& getMap() const { return contig_map_; }

private:

  VCFContigMap contig_map_;

  bool getCreateContig(const ContigId_t& contig_id, std::shared_ptr<ContigMultiMap>& contig_ptr);

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds multi ploid variants until they can be phased.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


using VCFGenomeMap = std::map<GenomeId_t, std::shared_ptr<VCFGenome>>;
class VCFPopulation {

public:

  explicit VCFPopulation() = default;
  VCFPopulation(const VCFPopulation&) = default;
  virtual ~VCFPopulation() = default;

  // Create the genome variant if it does not exist.
  bool getCreateGenome(const GenomeId_t& genome_id, std::shared_ptr<VCFGenome>& genome);

  size_t variantCount() const;

  const VCFGenomeMap& getMap() const { return genome_map_; }


private:

  VCFGenomeMap genome_map_;



};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object accepts multi ploid variants from the VCF parser and phases them.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


class GenomePhasing {

public:

  explicit GenomePhasing() = default;
  ~GenomePhasing() = default;

  static bool haploidPhasing(const VCFPopulation& vcf_population,
                             std::shared_ptr<const GenomeDatabase> genome_db,
                             std::shared_ptr<PopulationVariant> haploid_population);

private:



};












}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_KGL_VARIANT_FACTORY_VCF_PHASING_H
