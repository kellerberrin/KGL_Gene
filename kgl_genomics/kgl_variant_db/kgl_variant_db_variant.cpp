//
// Created by kellerberrin on 24/05/23.
//

#include "kgl_variant_db_variant.h"

#include <ranges>

namespace kgl = kellerberrin::genome;


void kgl::VariantDBVariant::createVariantDB(const std::shared_ptr<const PopulationDB>& population_ptr) {

  // Get the unique variants for this population.
  auto unique_variant_map = population_ptr->uniqueVariants();

  // Create the variant index.
  for (auto const& [index, hgvs_pair] : std::views::enumerate(unique_variant_map)) {

    auto const& [hgvs, variant_ptr] = hgvs_pair;
    auto [insert_iter, result] = variant_index_.try_emplace(hgvs, std::pair{variant_ptr, index});
    if (not result) {

      ExecEnv::log().error("VariantDBVariant::createVariantDB; Unable to insert variant: {}, unexpected duplicate", hgvs);

    }

  }
  const size_t variant_size = variant_index_.size();

  // Create the Genome index and genome data vector.
  for (auto const& [index, genome_pair] : std::views::enumerate(population_ptr->getMap())) {

    auto const& [genome_id, genome_ptr] = genome_pair;
    auto [insert_iter, result] = genome_index_.try_emplace(genome_id, index);
    if (not result) {

      ExecEnv::log().error("VariantDBVariant::createVariantDB; Unable to insert Genome: {}, unexpected duplicate", genome_id);
      continue;

    }

    genome_data_.emplace_back(genome_id, std::vector<uint8_t>(variant_size, 0));

  }

  // Multi-threaded update of the genome data vectors.
  // Thread safety: each concurrent task updates the (disjoint) variant vector of a single genome.
  population_ptr->processAll_MT([&](std::shared_ptr<const GenomeDB> genome_ptr, const std::shared_ptr<const Variant>& variant_ptr) -> bool {

    if (not genome_ptr or not variant_ptr) {

      ExecEnv::log().critical("VariantDBVariant::createVariantDB; Bad GenomeDB or Bad Variant pointer");

    }

    auto find_variant_iter = variant_index_.find(variant_ptr->HGVS());
    if (find_variant_iter == variant_index_.end()) {

      ExecEnv::log().error("VariantDBVariant::createVariantDB; Genome: {}, Variant: {} not found in variant index", genome_ptr->genomeId(), variant_ptr->HGVS());
      return false;

    }

    auto const& [var_hgvs, var_pair] = *find_variant_iter;
    auto const& [var_ptr, var_index] = var_pair;

    auto find_genome_iter = genome_index_.find(genome_ptr->genomeId());
    if (find_genome_iter == genome_index_.end()) {

      ExecEnv::log().error("VariantDBVariant::createVariantDB; Genome: {} not found in genome index", genome_ptr->genomeId());
      return false;

    }

    auto const& [genome, genome_index] = *find_genome_iter;

    // Bounds check before use (the genome data vector is structurally frozen before multi-threading commences).
    if (genome_index >= genome_data_.size()) {

      ExecEnv::log().critical("VariantDBVariant::createVariantDB; Genome: {}, Genome Index: {} exceeds genome data size: {}",
                              genome, genome_index, genome_data_.size());

    }

    auto& [vec_genome, variant_array] = genome_data_[genome_index];

    if (var_index >= variant_array.size()) {

      ExecEnv::log().critical("VariantDBVariant::createVariantDB; Genome: {}, Variant Index: {} exceeds variant vector size: {}",
                              genome, var_index, variant_array.size());

    }

    ++(variant_array[var_index]);

    return true;

  });

}


// Classify an allele count (0, 1, 2) into the allele summary. Logs a warning if the allele count is non-diploid.
void kgl::VariantDBVariant::classifyAllele(AlleleSummmary& allele_summary, uint8_t allele_type, const std::string& context) const {

  switch (allele_type) {

    case 0:
      ++allele_summary.referenceHomozygous_;
      break;

    case 1:
      ++allele_summary.minorHeterozygous_;
      break;

    case 2:
      ++allele_summary.minorHomozygous_;
      break;

    default:
      ExecEnv::log().warn("VariantDBVariant; {}, has non-diploid allele count: {}", context, allele_type);
      break;

  }

}


kgl::AlleleSummmary kgl::VariantDBVariant::summaryByVariant(const std::shared_ptr<const Variant>& variant) const {

  AlleleSummmary allele_summary;

  // Retrieve the variant index
  auto find_iter = variant_index_.find(variant->HGVS());
  if (find_iter == variant_index_.end()) {

    ExecEnv::log().error("VariantDBVariant::summaryByVariant; Unable to find variant: {}" , variant->HGVS());
    return allele_summary;

  }
  auto const& [hgvs, variant_pair] = *find_iter;
  auto const& [variant_ptr, variant_index] = variant_pair;

  // Loop through the Genomes.
  for (auto const& [genome, variant_vector] : genome_data_) {

    classifyAllele(allele_summary, variant_vector[variant_index], std::format("summaryByVariant; Genome: {}, Variant Index: {}", genome, variant_index));

  }

  // The sum of allele types should equal the number of Genomes.
  size_t allele_sum = allele_summary.referenceHomozygous_ + allele_summary.minorHomozygous_ + allele_summary.minorHeterozygous_;
  if (genome_data_.size() != allele_sum) {

    ExecEnv::log().warn("VariantDBVariant::summaryByVariant; Genome Vector Size: {}, Does not equal the sum of allele types: {}",
                        genome_data_.size() , allele_sum);

  }

  return allele_summary;

}


kgl::AlleleSummmary kgl::VariantDBVariant::summaryByGenome(const GenomeId_t& genome) const {

  AlleleSummmary allele_summary;

  // Retrieve the genome index
  auto find_iter = genome_index_.find(genome);
  if (find_iter == genome_index_.end()) {

    ExecEnv::log().error("VariantDBVariant::summaryByGenome; Unable to find genome: {}" , genome);
    return allele_summary;

  }
  auto const& [genome_id, genome_index] = *find_iter;
  auto const& [genome_, allele_vector] = genome_data_[genome_index];

  // Loop through the Variants.
  for (auto const& allele_type : allele_vector) {

    classifyAllele(allele_summary, allele_type, std::format("summaryByGenome; Genome: {}, Genome Index: {}", genome, genome_index));

  }

  size_t allele_sum = allele_summary.referenceHomozygous_ + allele_summary.minorHomozygous_ + allele_summary.minorHeterozygous_;
  if (variant_index_.size() != allele_sum) {

    ExecEnv::log().warn("VariantDBVariant::summaryByGenome; Variant Index Size: {}, Does not equal the sum of allele types: {}",
                        variant_index_.size() , allele_sum);

  }

  return allele_summary;

}


kgl::AlleleSummmary kgl::VariantDBVariant::populationSummary() const {

  AlleleSummmary allele_summary;

  // Retrieve through Genomes.
  for (auto const& [genome_id, allele_vector] : genome_data_) {

    // Loop through the Variants.
    for (auto const& allele_type : allele_vector) {

      classifyAllele(allele_summary, allele_type, std::format("populationSummary; Genome Id: {}", genome_id));

    }

  }

  size_t allele_sum = allele_summary.referenceHomozygous_ + allele_summary.minorHomozygous_ + allele_summary.minorHeterozygous_;
  size_t allele_type_count = variant_index_.size() * genome_data_.size();
  if (allele_type_count != allele_sum) {

    ExecEnv::log().warn("VariantDBVariant::populationSummary; Variants * Genomes: {}, Does not equal the sum of allele types: {}",
                        allele_type_count , allele_sum);

  }

  return allele_summary;

}