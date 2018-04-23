//
// Created by kellerberrin on 23/04/18.
//

#ifndef KGL_VARIANT_DB_UNPHASED_H
#define KGL_VARIANT_DB_UNPHASED_H



#include "kgl_utility.h"
#include "kgl_variant.h"


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

  // true if vector.size() < 2
  static bool isHomozygous(const std::vector<std::shared_ptr<Variant>>& variant_vector);

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

  // Generate phasing statitics.
  bool genomePhasingStats(const GenomeId_t& genome_id,
                          bool snp_only,
                          size_t& heterozygous,
                          size_t& homozygous,
                          size_t& singleheterozygous) const;

  bool getUnphasedVariants(const GenomeId_t& genome_id,
                           const ContigId_t& contig_id,
                           ContigOffset_t offset,
                           std::vector<std::shared_ptr<Variant>>& variant_vector) const;

private:

  UnphasedGenomeMap genome_map_;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_KGL_VARIANT_DB_UNPHASED_H
