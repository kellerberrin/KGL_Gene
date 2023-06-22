//
// Created by kellerberrin on 9/04/23.
//

#ifndef KGL_ANALYSIS_PFEMP_VARIANT_H
#define KGL_ANALYSIS_PFEMP_VARIANT_H


#include "kgl_genome_collection.h"
#include "kgl_variant_db_population.h"



namespace kellerberrin::genome {   //  organization::project level namespace



// A results class that holds the results of gene variant analysis for all genomes.
class GenomeGeneVariantAnalysis {

  using VariantGeneMap = std::map<std::string, std::pair<std::shared_ptr<const GeneFeature>, std::shared_ptr<ContigDB>>>;

public:

  GenomeGeneVariantAnalysis() {

    gene_population_ptr = std::make_shared<PopulationDB>("GenePopulation", DataSourceEnum::Falciparum);

  }
  ~GenomeGeneVariantAnalysis() = default;

  void writeGeneResults(const std::string &file_name);
  void setGeneVector(const GeneVector& gene_vector);
  void getGeneVariants(const std::shared_ptr<const PopulationDB> &population_ptr);

private:

  GeneVector gene_vector_;   // Gene vector of interest.
  std::shared_ptr<PopulationDB> gene_population_ptr; // Each contig per genome is a gene in the gene vector

  constexpr static const char CSV_DELIMITER_ = ',';


};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The genome statistics for each variant.
struct GenomeCount {

  std::shared_ptr<const Variant> variant_ptr_;
  std::set<GenomeId_t> genome_set_;
  std::set<GenomeId_t> homozygous_set_;

};
using GenomeCountMap = std::map<std::string, std::shared_ptr<GenomeCount>>;
using GenomeCountSorted = std::multimap<size_t, std::shared_ptr<const GenomeCount>, std::greater<>>;

class GeneGenomeAnalysis {


public:


  // Initialized with the unique variants in the coding region of a Gene.
  GeneGenomeAnalysis(std::shared_ptr<const GeneFeature> gene_ptr,
                     const std::shared_ptr<const ContigDB>& gene_unique_variants);
  ~GeneGenomeAnalysis() = default;

  void analyzeGenePopulation(const std::shared_ptr<const PopulationDB>& gene_population_ptr);

  [[nodiscard]] const std::vector<GenomeId_t>& zeroVariants() const { return zero_variants_; }

  // Returns the GenomeCount objects in descending count order.
  [[nodiscard]] GenomeCountSorted getCountSorted() const;
  // Returns the variants present in only one genome (singletons)
  [[nodiscard]] std::vector<std::shared_ptr<const GenomeCount>> getSingletonVariants() const;
  // Returns the genomes containing variants unique to that genome (singleton genomes).
  // Note that this is not the same as getSingletonVariants() since we can have multiple singleton variants
  // per singleton genome. In general, getSingletonGenomes() <= getSingletonVariants().
  [[nodiscard]] std::set<kgl::GenomeId_t> getSingletonGenomes() const;
  // Returns the top n variant counts.
  template<size_t N>
  [[nodiscard]] std::array<size_t, N> getTopCounts() const {

    std::array<size_t, N> top_count_variants{0};
    size_t index = 0;
    for (auto const& [count, genome_count] : getCountSorted()) {

      if (index >= N) {

        break;

      }

      top_count_variants[index] = count;

      ++index;

    }

    return top_count_variants;

  }

private:

  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<GenomeCountMap> gene_genome_analysis_ptr_;
  std::vector<GenomeId_t> zero_variants_;

  bool addVariant(const std::shared_ptr<const Variant>& variant_ptr);

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct GenotypeCount {

  std::shared_ptr<const ContigDB> genotype_;
  std::vector<GenomeId_t> genomes_;
  size_t genotype_hash_;

};

using GenotypeMap = std::map<size_t, std::shared_ptr<GenotypeCount>>;
using GenotypeCountSorted = std::multimap<size_t, std::shared_ptr<const GenotypeCount>, std::greater<>>;

class GenotypeAnalysis {

public:

  explicit GenotypeAnalysis(std::shared_ptr<const GeneFeature> gene_ptr) : gene_ptr_(std::move(gene_ptr)) {}
  ~GenotypeAnalysis() = default;

  void analyzeGenePopulation(const std::shared_ptr<const PopulationDB>& gene_population_ptr);
  [[nodiscard]] const GenotypeMap& getGenoTypeMap() const { return genotype_map_; }
  // Returns the GenomeCount objects in descending count order.
  [[nodiscard]] GenotypeCountSorted getCountSorted() const;

  // Returns the top n genotypes.
  template<size_t N>
  [[nodiscard]] std::array<size_t, N> getTopCounts() const {

    std::array<size_t, N> top_count_genotypes{0};
    size_t index = 0;
    for (auto const& [count, genotype_ptr] : getCountSorted()) {

      if (index >= N) {

        break;

      }

      top_count_genotypes[index] = count;

      ++index;

    }

    return top_count_genotypes;

  }

  // Quantitative Variation Statistics for category (genotype) data.
  [[nodiscard]] double Gini() const;   // The Gini coefficient.
  [[nodiscard]] double HRel() const;   // A modified Shannon entropy measure for different category (genotype) counts.
  [[nodiscard]] double IQV() const;    // Index of qualitative variation.
  [[nodiscard]] double IQVAdj() const { return (2.0-IQV()); }    // Adjusted Index of qualitative variation.


private:

  std::shared_ptr<const GeneFeature> gene_ptr_;
  GenotypeMap genotype_map_;
  std::vector<GenomeId_t> zero_variants_;

  size_t genotypeHash(const std::shared_ptr<const ContigDB>& genotype);


};



}; // Namespace.


#endif // KGL_ANALYSIS_PFEMP_VARIANT_H
