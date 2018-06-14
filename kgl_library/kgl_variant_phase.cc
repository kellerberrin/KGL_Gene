//
// Created by kellerberrin on 23/04/18.
//

#include "kgl_variant_phase.h"
#include "kgl_variant_evidence.h"
#include "kgl_variant_single.h"
#include <fstream>


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Not Thread safe.
// This object accepts unphased variants and trivally phases them as haploid.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


// Just copy into a population object.
bool kgl::GenomePhasing::haploidPhasing(std::shared_ptr<const UnphasedPopulation> unphased_population_ptr,
                                        std::shared_ptr<const GenomeDatabase> genome_db,
                                        std::shared_ptr<PhasedPopulation> haploid_population)  {

  bool result = true;

  for (auto genome : unphased_population_ptr->getMap()) {

    // Create the GenomeVariant object.
    std::shared_ptr<GenomeVariant> genome_variant = GenomeVariant::emptyGenomeVariant(genome.first, GenomeVariant::HAPLOID_GENOME, genome_db);

    // Add all the contig.
    for (auto contig : genome.second->getMap()) {

      // Add all offsets.
      for (auto offset : contig.second->getMap()) {

        // add all the variants
        for (auto variant : offset.second) {

          std::shared_ptr<Variant> mutable_variant = std::const_pointer_cast<Variant>(variant);
          mutable_variant->phaseId(ContigVariant::HAPLOID_HOMOLOGOUS_INDEX);   // Assign to the first (and only) homologous contig.
          genome_variant->addVariant(variant);

        } // All variants.

      } // All offsets

    } // All contigs

    if (not haploid_population->addGenomeVariant(genome_variant)) {

      ExecEnv::log().error("Unable to add genome: {} to haploid population", genome_variant->genomeId());
      result = false;

    }

  } // All genomes

  ExecEnv::log().info("Haploid Phasing, pre_phase variants: {}, genomes: {}, resultant population variants: {}, genomes: {}",
                      unphased_population_ptr->variantCount(), unphased_population_ptr->getMap().size(),
                      haploid_population->variantCount(), haploid_population->getMap().size());

  return result;

}




////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These objects accept unphased variants and writes out ref/alt read statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::VariantClassifier::VariantClassifier(std::shared_ptr<const UnphasedPopulation> vcf_population_ptr) {

  // Get a list of all genomes.
  genomes_ = vcf_population_ptr->genomeList();

  // for all variants index by contig/offset
  for (auto genome : vcf_population_ptr->getMap()) {

    for(auto contig : genome.second->getMap()) {

      for(auto offset : contig.second->getMap()) {

        for(auto variant : offset.second) {

          ContigOffsetVariants ordered_variant(variant);
          auto result = variant_map_.find(ordered_variant);

          if (result != variant_map_.end()) {


            // We ignore duplicate insert errors, since this is just the variant caller calling
            // 2 identical variants under the ploidy = 2 setting. Both variants have identical ref/alt ratios.
            // with a small (generally zero) ref count and a much larger alt count.
            result->second.insert(GenomeOffsetMap::value_type(variant->genomeId(), variant));

          } else {

            GenomeOffsetMap genome_map;

            auto insert_result = variant_map_.insert(VariantMap::value_type(ordered_variant, genome_map));

            if (not insert_result.second) {

              ExecEnv::log().error("Attempt to insert duplicate variant: {}",
                                   variant->output(' ',VariantOutputIndex::START_0_BASED, false));

            } else {

              // We ignore duplicate insert errors, since this is just the variant caller calling
              // 2 identical variants under the ploidy = 2 setting. Both variants have identical ref/alt ratios.
              // with a small (generally zero) ref count and a much larger alt count.
              insert_result.first->second.insert(GenomeOffsetMap::value_type(variant->genomeId(), variant));

            }

          } // if found

        } // for variant

      } // for offset

    } // for contig

  } // for genome

}


bool kgl::VariantClassifier::writeVariants(char delimiter,
                                           const std::string& file_name,
                                           size_t min_count,
                                           bool ref /* false is alt*/) const {

  std::ofstream variant_file(file_name);

  if (not variant_file.good()) {

    ExecEnv::log().error("Unable to open variant file: {}", file_name);
    return false;

  }

  // write the header.
  variant_file << "Contig" << delimiter << "Offset" << delimiter;

  for (auto genome : getGenomes()) {

    variant_file << genome << delimiter;

  }
  variant_file << std::endl;

  // for all variant offsets.
  for (auto variant_offset : getMap()) {

    variant_file << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset() << delimiter;

    for (auto genome : getGenomes()) {

      auto find_result = variant_offset.second.find(genome);

      if (find_result != variant_offset.second.end()) {

        std::shared_ptr<const CountEvidence> count_evidence_ptr = std::dynamic_pointer_cast<const CountEvidence>(find_result->second->evidence());

        if (count_evidence_ptr) {

          size_t total_count = count_evidence_ptr->refCount() + count_evidence_ptr->altCount();

          if (total_count >= min_count) {

            if (ref) {

              variant_file << count_evidence_ptr->refCount() << delimiter;

            } else {

              variant_file << count_evidence_ptr->altCount() << delimiter;

            }

          } else {

            variant_file << '0' << delimiter;

          }

        } else {

          ExecEnv::log().error("Variant without count evidence: {}", find_result->second->output(' ', VariantOutputIndex::START_0_BASED, false));

        }

      } else { // genome not found

        variant_file << '0' << delimiter;

      }

    }

    variant_file << std::endl;

  }

  return true;

}


