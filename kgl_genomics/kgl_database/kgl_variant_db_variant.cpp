//
// Created by kellerberrin on 24/05/23.
//

#include "kgl_variant_db_variant.h"


namespace kgl = kellerberrin::genome;


void kgl::VariantDBVariant::createVariantDB(const std::shared_ptr<const PopulationDB>& population_ptr) {

  // Create the variant index.
  auto variant_index_ptr = std::make_unique<VariantDBVariantIndex>();
  auto unique_variant_map = population_ptr->uniqueVariants();
  size_t index{0};
  for (auto &[hgvs, variant_ptr]: unique_variant_map) {

    auto [insert_iter, result] = variant_index_ptr->try_emplace(hgvs, std::pair<std::shared_ptr<const Variant>, size_t>{variant_ptr, index});
    if (not result) {

      ExecEnv::log().error("VariantDBVariant::createVariantDB; Unable to insert variant: {}, unexpected duplicate", hgvs);
      continue;

    }

    ++index;

  }

  size_t variant_size = variant_index_.size();

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

      ++variant_array[var_index];

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


