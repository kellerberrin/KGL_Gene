//
// Created by kellerberrin on 24/05/23.
//

#ifndef KGL_VARIANT_DB_VARIANT_H
#define KGL_VARIANT_DB_VARIANT_H

#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// VariantDBVariant; Indexed by Variant and Genome. Genomes containing the particular variant (Hom/Het counted) are listed.
// In addition (importantly), each Genome the does NOT contain the population variant is listed. In other words, Genomes that carry
// the homozygous reference allele instead of the variant are also listed for the particular variant.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Auxiliary class to return allele statistics by Genome, Variant or Population.
// If queried by Variant, the sum of the summary terms should equal the number of Genomes.
// If queried by Genome, the sum of the summary terms should equal the number of Variants.
// If a population summary, the sum of the summary terms should equal the number of Variants x the number of Genomes.
struct AlleleSummmary {

  size_t referenceHomozygous_{0};    // (A,A)
  size_t minorHeterozygous_{0};      // (a, A) [or rarely can be (a, b)]
  size_t minorHomozygous_{0};        // (a, a)

  // Convenience operator
  void operator+=(const AlleleSummmary& rhs) {

    minorHomozygous_ += rhs.minorHomozygous_;
    referenceHomozygous_ += rhs.referenceHomozygous_;
    minorHeterozygous_ += rhs.minorHeterozygous_;

  }

};

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

  [[nodiscard]] AlleleSummmary summaryByVariant(const std::shared_ptr<const Variant>& variant) const;
  [[nodiscard]] AlleleSummmary summaryByGenome(const GenomeId_t& genome) const;
  [[nodiscard]] AlleleSummmary populationSummary() const;

private:

  VariantDBGenomeIndex genome_index_;
  VariantDBVariantIndex variant_index_;
  VariantDBGenomeData genome_data_;

  void createVariantDB(const std::shared_ptr<const PopulationDB>& population_ptr);

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Given a contig and offset address (or just a variant_ptr).
// This object returns the feature record (gene feature) of the nearest coding sequence.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////

class FindCodingInterval; // Pimpl implementation for determining variant membership of gene coding intervals.

class FindFeature {

public:

  explicit FindFeature(const std::shared_ptr<const GenomeReference>& reference_ptr);
  ~FindFeature();

private:

  std::unique_ptr<const FindCodingInterval> pimpl_coding_interval_ptr_;

};



} // end namespace

#endif //KGL_VARIANT_DB_VARIANT_H
