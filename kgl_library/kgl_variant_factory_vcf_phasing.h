//
// Created by kellerberrin on 2/04/18.
//

#ifndef KGL_VARIANT_FACTORY_VCF_PHASING_H
#define KGL_VARIANT_FACTORY_VCF_PHASING_H


#include "kgl_utility.h"
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants until they can be phased.
// This object hold variants for each contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using UnphasedOffsetMap = std::map<ContigOffset_t, std::vector<std::shared_ptr<Variant>>>;
class UnphasedContig {

public:

  explicit UnphasedContig(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  UnphasedContig(const UnphasedContig&) = default;
  virtual ~UnphasedContig() = default;

  const ContigId_t& contigId() const { return contig_id_; }

  bool addVariant(std::shared_ptr<Variant> variant);

  size_t variantCount() const;

  const UnphasedOffsetMap& getMap() const { return contig_offset_map_; }

private:

  ContigId_t contig_id_;
  UnphasedOffsetMap contig_offset_map_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants until they can be phased.
// This object hold variants for each genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using UnphasedContigMap = std::map<ContigId_t, std::shared_ptr<UnphasedContig>>;
class UnphasedGenome {

public:

  explicit UnphasedGenome(const GenomeId_t& genome_id) : genome_id_(genome_id) {}
  UnphasedGenome(const UnphasedGenome&) = default;
  virtual ~UnphasedGenome() = default;

  size_t variantCount() const;

  bool addVariant(std::shared_ptr<Variant> variant);

  const GenomeId_t& genomeId() const { return genome_id_; }

  const UnphasedContigMap& getMap() const { return contig_map_; }

private:

  UnphasedContigMap contig_map_;
  GenomeId_t genome_id_;

  bool getCreateContig(const ContigId_t& contig_id, std::shared_ptr<UnphasedContig>& contig_ptr);

};


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
  UnphasedPopulation(const UnphasedPopulation&) = default;
  virtual ~UnphasedPopulation() = default;

  // Create the genome variant if it does not exist.
  bool getCreateGenome(const GenomeId_t& genome_id, std::shared_ptr<UnphasedGenome>& genome);

  size_t variantCount() const;

  const UnphasedGenomeMap& getMap() const { return genome_map_; }


private:

  UnphasedGenomeMap genome_map_;



};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object accepts unphased variants from the VCF parser and phases them.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


class GenomePhasing {

public:

  explicit GenomePhasing() = default;
  ~GenomePhasing() = default;

  static bool haploidPhasing(std::shared_ptr<const UnphasedPopulation> vcf_population_ptr,
                             std::shared_ptr<const GenomeDatabase> genome_db,
                             std::shared_ptr<PhasedPopulation> haploid_population);

private:



};





}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_KGL_VARIANT_FACTORY_VCF_PHASING_H
