//
// Created by kellerberrin on 23/04/18.
//

#ifndef KGL_VARIANT_PHASE_H
#define KGL_VARIANT_PHASE_H



#include "kgl_utility.h"
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



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


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These objects accept unphased variants and writes out ref/alt read statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigOffsetVariants {

public:

  explicit ContigOffsetVariants(std::shared_ptr<const Variant> variant) : variant_(variant) {}
  explicit ContigOffsetVariants(const ContigOffsetVariants& copy) = default;
  ~ContigOffsetVariants() = default;

  std::shared_ptr<const Variant> variant() const { return variant_; }

  bool operator<(const ContigOffsetVariants& cmp) const { return variant_->lessThan(*cmp.variant()); }

private:

  std::shared_ptr<const Variant> variant_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Stores all the variants in contig/offset order.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeOffsetMap = std::map<GenomeId_t, std::shared_ptr<const Variant>>;
using VariantMap = std::map<ContigOffsetVariants, GenomeOffsetMap>;
class VariantClassifier {

public:

  explicit VariantClassifier(std::shared_ptr<const UnphasedPopulation> vcf_population_ptr);
  ~VariantClassifier() = default;

  const VariantMap& getMap() const { return variant_map_; }
  const std::vector<GenomeId_t>& getGenomes() const { return genomes_; }

  bool writeVariants(char delimiter,
                     const std::string& file_name,
                     size_t min_count,
                     bool ref /* false is alt*/) const;

private:

  VariantMap variant_map_;
  std::vector<GenomeId_t> genomes_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_PHASE_H
