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
// This object holds variants for each contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////



using UnphasedVariantCount = std::pair<std::shared_ptr<const Variant>, size_t>;
using UnphasedVectorVariantCount = std::vector<UnphasedVariantCount>;
using UnphasedOffsetMap = std::map<ContigOffset_t, UnphasedVectorVariantCount>;
class UnphasedContig {

public:

  explicit UnphasedContig(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  UnphasedContig(const UnphasedContig&) = delete;
  virtual ~UnphasedContig() = default;

  UnphasedContig& operator=(const UnphasedContig&) = delete; // Use deep copy.

  const ContigId_t& contigId() const { return contig_id_; }

  bool addVariant(std::shared_ptr<const Variant> variant);

  size_t variantCount() const;

  const UnphasedOffsetMap& getMap() const { return contig_offset_map_; }

  std::shared_ptr<UnphasedContig> filterVariants(const VariantFilter& filter) const;

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
  UnphasedGenome(const UnphasedGenome&) = delete;
  virtual ~UnphasedGenome() = default;

  UnphasedGenome& operator=(const UnphasedGenome&) = delete; // Use deep copy.

  size_t variantCount() const;

  bool addVariant(std::shared_ptr<const Variant> variant);

  const GenomeId_t& genomeId() const { return genome_id_; }

  std::shared_ptr<UnphasedGenome> filterVariants(const VariantFilter& filter) const;

  const UnphasedContigMap& getMap() const { return contig_map_; }

  bool getCreateContig(const ContigId_t& contig_id, std::shared_ptr<UnphasedContig>& contig_ptr);

private:

  UnphasedContigMap contig_map_;
  GenomeId_t genome_id_;

  bool addContig(std::shared_ptr<UnphasedContig> contig_ptr);

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_KGL_VARIANT_DB_UNPHASED_H
