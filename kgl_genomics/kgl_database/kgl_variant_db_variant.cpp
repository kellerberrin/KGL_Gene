//
// Created by kellerberrin on 24/05/23.
//

#include "kgl_variant_db_variant.h"


namespace kgl = kellerberrin::genome;






void kgl::VariantDBVariant::createVariantDB(const std::shared_ptr<const PopulationDB>& population_ptr) {

  VariantGenomeMap genome_map;
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    genome_map[genome_id] = 0;

  }

  // Add a local class to process each variant.
  class InsertVariants {

  public:

    InsertVariants(const VariantGenomeMap& genome_map) : genome_map_(genome_map) {}

    bool addVariants(const std::shared_ptr<const Variant>& variant) {

      if (not variant_map_.contains(variant->HGVS())) {


        auto [insert_iter, result] = variant_map_.try_emplace(variant->HGVS(), VariantDBRecord(variant, genome_map_));
        if (not result) {

          ExecEnv::log().error("VariantDBVariant::createVariantDB; Unable to insert variant: {}, unexpected duplicate", variant->HGVS());

        }

      }

      return true;

    }

    bool countVariants(const std::shared_ptr<const Variant>& variant) {

      auto variant_hgvs = variant->HGVS();

      if (not variant_map_.contains(variant_hgvs)) {

        ExecEnv::log().error("VariantDBVariant::createVariantDB; variant: {} not found in map.", variant_hgvs);
        return false;

      }

      auto& [hgvs, variant_record] = *variant_map_.find(variant_hgvs);

      if (genome_.empty()) {

        ExecEnv::log().error("VariantDBVariant::createVariantDB; genome not initialized for variant: {}", variant_hgvs);
        return false;

      }

      auto genome_iter = variant_record.getMap().find(genome_);
      if (genome_iter == variant_record.getMap().end()) {

        ExecEnv::log().error("VariantDBVariant::createVariantDB; genome: {} not found for variant: {}", genome_, variant_hgvs);
        return false;

      }

      auto& [genome_id, variant_count] = *genome_iter;
      ++variant_count;

      return true;

    }

    VariantGenomeMap genome_map_;
    VariantVariantMap variant_map_;
    GenomeId_t genome_;

  };

  InsertVariants insert_variants(genome_map);
  // For all variants populate with zero counted genomes.
  population_ptr->processAll(insert_variants, &InsertVariants::addVariants);

  // Count the all variant/genomes.
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    insert_variants.genome_ = genome_id; // Must setup the genome id within the local class, see the InsertVariants::countVariants logic.
    genome_ptr->processAll(insert_variants, &InsertVariants::countVariants);

  }

  // Move the local class map to the VariantDBVariant map.
  variant_map_ = std::move(insert_variants.variant_map_);

}
