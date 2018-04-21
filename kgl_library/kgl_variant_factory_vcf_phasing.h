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
// An internal parser variant object that holds multi ploid variants until they can be phased.
// This object hold variants for each contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using ContigOffsetMap = std::map<ContigOffset_t, std::vector<std::shared_ptr<Variant>>>;
class VCFContig {

public:

  explicit VCFContig(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  VCFContig(const VCFContig&) = default;
  virtual ~VCFContig() = default;

  const ContigId_t& contigId() const { return contig_id_; }

  bool addVariant(std::shared_ptr<Variant> variant);

  size_t variantCount() const;

  const ContigOffsetMap& getMap() const { return contig_offset_map_; }

private:

  ContigId_t contig_id_;
  ContigOffsetMap contig_offset_map_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds multi ploid variants until they can be phased.
// This object hold variants for each genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using VCFContigMap = std::map<ContigId_t, std::shared_ptr<VCFContig>>;
class VCFGenome {

public:

  explicit VCFGenome(const GenomeId_t& genome_id) : genome_id_(genome_id) {}
  VCFGenome(const VCFGenome&) = default;
  virtual ~VCFGenome() = default;

  size_t variantCount() const;

  bool addVariant(std::shared_ptr<Variant> variant);

  const GenomeId_t& genomeId() const { return genome_id_; }

  const VCFContigMap& getMap() const { return contig_map_; }

private:

  VCFContigMap contig_map_;
  GenomeId_t genome_id_;

  bool getCreateContig(const ContigId_t& contig_id, std::shared_ptr<VCFContig>& contig_ptr);

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds multi ploid variants until they can be phased.
// This object hold variants for a population.
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

  static bool haploidPhasing(std::shared_ptr<const VCFPopulation> vcf_population_ptr,
                             std::shared_ptr<const GenomeDatabase> genome_db,
                             std::shared_ptr<PopulationVariant> haploid_population);

private:



};





}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_KGL_VARIANT_FACTORY_VCF_PHASING_H
