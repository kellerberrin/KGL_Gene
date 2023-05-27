//
// Created by kellerberrin on 24/05/23.
//

#ifndef KGL_VARIANT_DB_VARIANT_H
#define KGL_VARIANT_DB_VARIANT_H

#include "kgl_variant_db_population.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These objects are NOT part of the PopulationDB structure. They exist to reorganise a population of
// Variants and Genomes into a data structure that can be conveniently queried such that:
//
// 1. VariantDBVariant; Indexed by Variant. Genomes containing a particular population variant (Hom/Het counted) are listed.
// In addition (importantly), each Genome the does NOT contain the population variant are listed. In other words, Genomes that carry
// the homozygous reference allele instead of the population variant are also listed for the particular variant.
//
// 2. VariantDBGenome; Indexed by Genome. All variants (Hom/Het counted) contained in a Genome are listed.
// In addition (importantly), all variants in the reference population that are not present in the Genome are also listed.
// In other words, where the Genome contains homozygous reference alleles at the Contig/Offset location where
// there are variants present in the reference population.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// VariantDBVariant; Indexed by Variant. Genomes containing the particular variant (Hom/Het counted) are listed.
// In addition (importantly), each Genome the does NOT contain the population variant is listed. In other words, Genomes that carry
// the homozygous reference allele instead of the variant are also listed for the particular variant.//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Indexed by variant HGVS signature, .second of the pair is the offset within the data vector.
using VariantDBVariantIndex = std::map<std::string, std::pair<std::shared_ptr<const Variant>, size_t>>;

// Allows the VariantDBGenomeData structure to be indexed directly using GenomeId_t.
using VariantDBGenomeIndex = std::map<GenomeId_t , size_t>;

// Each element in the data vector is valued either as; no allele (0), heterozygous(1) or homozygous(2).
// The offset within the vector corresponds to the offset in the VariantDBVariantIndex map (the variant size_t offset).
using VariantDBGenomeData = std::vector<std::pair<GenomeId_t , std::vector<uint8_t>>>;

class VariantDBVariant {

public:

  explicit VariantDBVariant(const std::shared_ptr<const PopulationDB>& population_ptr) { createVariantDB(population_ptr); }
  ~VariantDBVariant() = default;

  [[nodiscard]] const VariantDBGenomeIndex& genomeMap() const { return genome_index_; }
  [[nodiscard]] const VariantDBVariantIndex& variantMap() const { return variant_index_; }
  [[nodiscard]] const VariantDBGenomeData& genomeData() const { return genome_data_; }

private:

  VariantDBGenomeIndex genome_index_;
  VariantDBVariantIndex variant_index_;
  VariantDBGenomeData genome_data_;

  void createVariantDB(const std::shared_ptr<const PopulationDB>& population_ptr);

};




} // end namespace

#endif //KGL_VARIANT_DB_VARIANT_H
