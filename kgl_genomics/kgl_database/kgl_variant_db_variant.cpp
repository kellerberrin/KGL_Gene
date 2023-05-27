//
// Created by kellerberrin on 24/05/23.
//

#include "kgl_variant_db_variant.h"


namespace kgl = kellerberrin::genome;


void kgl::VariantDBVariant::createVariantDB(const std::shared_ptr<const PopulationDB>& population_ptr) {

  // Create the variant index.
  auto unique_variant_map = population_ptr->uniqueVariants();
  size_t index{0};
  for (auto &[hgvs, variant_ptr]: unique_variant_map) {

    auto [insert_iter, result] = variant_index_.try_emplace(hgvs, std::pair<std::shared_ptr<const Variant>, size_t>{variant_ptr, index});
    if (not result) {

      ExecEnv::log().error("VariantDBVariant::createVariantDB; Unable to insert variant: {}, unexpected duplicate", hgvs);
      continue;

    }

    ++index;

  }

  size_t variant_size = variant_index_.size();

  // Create the Genome index and std>::shared_ptr data vector
  std::shared_ptr<VariantDBGenomeData> genome_data_ptr;
  index = 0;
  for (auto &[genome_id, genome_ptr]: population_ptr->getMap()) {

    auto [insert_iter, result] = genome_index_.try_emplace(genome_id, index);
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

    LocalGenomeVariant(std::shared_ptr<VariantDBGenomeData> genome_data_ptr) : genome_data_ptr_(std::move(genome_data_ptr)) {}

    bool processGenomeVariants(std::shared_ptr<const GenomeDB> genome_ptr, const std::shared_ptr<const Variant>& variant_ptr) {


      return true;

    }

    std::shared_ptr<VariantDBGenomeData> genome_data_ptr_;

  };

  LocalGenomeVariant genome_variant(genome_data_ptr);

  population_ptr->processAll_MT(genome_variant, &LocalGenomeVariant::processGenomeVariants);

  genome_data_ = std::move(*genome_variant.genome_data_ptr_);

}


