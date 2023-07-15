//
// Created by kellerberrin on 24/05/23.
//

#include "kgl_variant_db_variant.h"


namespace kgl = kellerberrin::genome;


void kgl::VariantDBVariant::createVariantDB(const std::shared_ptr<const PopulationDB>& population_ptr) {

  // Get the unique variants for this population.
  auto unique_variant_map = population_ptr->uniqueVariants();
  // Create the variant index.
  auto variant_index_ptr = std::make_unique<VariantDBVariantIndex>();
  size_t index{0};
  for (auto &[hgvs, variant_ptr]: unique_variant_map) {

    auto [insert_iter, result] = variant_index_ptr->try_emplace(hgvs, std::pair<std::shared_ptr<const Variant>, size_t>{variant_ptr, index});
    if (not result) {

      ExecEnv::log().error("VariantDBVariant::createVariantDB; Unable to insert variant: {}, unexpected duplicate", hgvs);
      continue;

    }

    ++index;

  }
  const size_t variant_size = variant_index_ptr->size();

  // Create the Genome index and std>::shared_ptr data vector
  auto genome_index_ptr = std::make_unique<VariantDBGenomeIndex>();
  auto genome_data_ptr = std::make_unique<VariantDBGenomeData>();
  index = 0;
  for (auto &[genome_id, genome_ptr]: population_ptr->getMap()) {

    auto [insert_iter, result] = genome_index_ptr->try_emplace(genome_id, index);
    if (not result) {

      ExecEnv::log().error("VariantDBVariant::createVariantDB; Unable to insert Genome: {}, unexpected duplicate", genome_id);
      continue;

    }

    genome_data_ptr->emplace_back(genome_id, std::vector<uint8_t>(variant_size, 0));

    ++index;

  }

  // Local class to update the genome data vectors.
  class LocalGenomeVariant {

  public:

    LocalGenomeVariant(std::unique_ptr<VariantDBVariantIndex> variant_index_ptr,
                       std::unique_ptr<VariantDBGenomeIndex> genome_index_ptr,
                       std::unique_ptr<VariantDBGenomeData> genome_data_ptr)
                       : variant_index_ptr_(std::move(variant_index_ptr)),
                       genome_index_ptr_(std::move(genome_index_ptr)),
                       genome_data_ptr_(std::move(genome_data_ptr)) {}

    bool processGenomeVariants(std::shared_ptr<const GenomeDB> genome_ptr, const std::shared_ptr<const Variant>& variant_ptr) {

      if (not genome_ptr or not variant_ptr) {

        ExecEnv::log().critical("LocalGenomeVariant::processGenomeVariants; Bad GenomeDB or Bad Variant pointer");

      }

      auto hgvs = variant_ptr->HGVS();
      auto find_variant_iter = variant_index_ptr_->find(hgvs);
      if (find_variant_iter == variant_index_ptr_->end()) {

        ExecEnv::log().error("VariantDBVariant::createVariantDB; Genome: {}, Variant: {} not found in variant index", genome_ptr->genomeId(), hgvs);
        return false;

      }

      auto const& [var_hgvs, var_pair] = *find_variant_iter;
      auto const& [var_ptr, var_index] = var_pair;

      auto find_genome_iter = genome_index_ptr_->find(genome_ptr->genomeId());
      if (find_genome_iter == genome_index_ptr_->end()) {

        ExecEnv::log().error("VariantDBVariant::createVariantDB; Genome: {} not found in genome index", genome_ptr->genomeId());
        return false;

      }

      auto const& [genome, genome_index] = *find_genome_iter;
      auto& [vec_genome, variant_array] = genome_data_ptr_->at(genome_index);

      if (var_index >= variant_array.size() or genome_index >= genome_data_ptr_->size()) {

        ExecEnv::log().critical("VariantDBVariant::createVariantDB; Genome: {}, Genome Index: {}, Variant Index: {} exceeds variant vector size: {}",
                                genome, genome_index, var_index, variant_array.size());

      }

      ++(variant_array[var_index]);

      return true;

    }

    std::unique_ptr<VariantDBVariantIndex> variant_index_ptr_;
    std::unique_ptr<VariantDBGenomeIndex> genome_index_ptr_;
    std::shared_ptr<VariantDBGenomeData> genome_data_ptr_;

  };


  LocalGenomeVariant genome_variant(std::move(variant_index_ptr), std::move(genome_index_ptr), std::move(genome_data_ptr));
  population_ptr->processAll_MT(genome_variant, &LocalGenomeVariant::processGenomeVariants);

  variant_index_ = std::move(*genome_variant.variant_index_ptr_);
  genome_index_ = std::move(*genome_variant.genome_index_ptr_);
  genome_data_ = std::move(*genome_variant.genome_data_ptr_);

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

    switch (variant_vector[variant_index]) {

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
        ExecEnv::log().warn("VariantDBVariant::summaryByVariant; Genome: {}, Variant Index: {} has non diploid allele count: {}",
                            genome, variant_index, variant_vector[variant_index]);
        break;

    }

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

  // Retrieve the variant index
  auto find_iter = genome_index_.find(genome);
  if (find_iter == genome_index_.end()) {

    ExecEnv::log().error("VariantDBVariant::summaryByGenome; Unable to find genome: {}" , genome);
    return allele_summary;

  }
  auto const& [genome_id, genome_index] = *find_iter;
  auto const& [genome_, allele_vector] = genome_data_[genome_index];

  // Loop through the Variants.
  for (auto const& allele_type : allele_vector) {

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
        ExecEnv::log().warn("VariantDBVariant::summaryByGenome; Genome: {}, Genome Index: {} has non-diploid allele count: {}",
                            genome, genome_index, allele_type);
        break;

    }

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
    for (auto const &allele_type: allele_vector) {

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
          ExecEnv::log().warn("VariantDBVariant::populationSummary; Genome Id: {}, has non-diploid allele count: {}", genome_id, allele_type);
          break;

      }

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
