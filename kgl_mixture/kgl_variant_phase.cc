//
// Created by kellerberrin on 23/04/18.
//

#include "kgl_variant_phase.h"
#include "kgl_read_phasing.h"
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

          std::shared_ptr<Variant> mutable_variant = std::const_pointer_cast<Variant>(variant.first);
          mutable_variant->updatePhaseId(ContigVariant::HAPLOID_HOMOLOGOUS_INDEX);   // Assign to the first (and only) homologous contig.
          genome_variant->addVariant(mutable_variant);

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



///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Not Thread safe.
// This object accepts a file with phasing information (Pf3k only) and phases the variants for a haploid
// organism. Important - only samples with a Mixture coefficient of 1 are phased. Samples with
// multiple strains are ignored.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::GenomePhasing::fileHaploidPhasing(const std::string& phase_file,
                                            size_t vcf_ploidy,
                                            std::shared_ptr<const UnphasedPopulation> unphased_population_ptr,
                                            std::shared_ptr<const GenomeDatabase> genome_db,
                                            std::shared_ptr<PhasedPopulation> haploid_population)  {

  bool result = true;

  GenomeMixtureStatistics mixture_statistics;

  if (not mixture_statistics.readMixtureStatistics(phase_file)) {

    ExecEnv::log().info("GenomePhasing::fileHaploidPhasing; cannot read phasing file: {}", phase_file);
    return false;

  }

  size_t statistics_count = 0;
  size_t mono_strain_count = 0;
  size_t heterozygous_count = 0;
  size_t accepted_heterozygous_count = 0;
  size_t homozygous_count = 0;
  for (auto genome : unphased_population_ptr->getMap()) {

    std::shared_ptr<GenomeVariant> genome_variant;
    if (not haploid_population->getGenomeVariant(genome.first, genome_variant)) {

      genome_variant = GenomeVariant::emptyGenomeVariant(genome.first, GenomeVariant::HAPLOID_GENOME, genome_db);
      if (not haploid_population->addGenomeVariant(genome_variant)) {

        ExecEnv::log().error("Haploid Phasing; Unable to add genome: {} to haploid population", genome_variant->genomeId());
        return false;

      }

    }

    // Create the GenomeVariant object.

    MixtureStatistics genome_statistics;
    if (not mixture_statistics.getMixtureStatistics(genome.first, genome_statistics)) {

      ExecEnv::log().info("GenomePhasing::fileHaploidPhasing; No phasing statistics available for genome: {}", genome.first);
      continue;

    }

    ++statistics_count;

    if (genome_statistics.first > 1) {

      ExecEnv::log().info("GenomePhasing::fileHaploidPhasing; Strains: {} mix for genome: {}", genome_statistics.first, genome.first);
      continue;

    }

    ++mono_strain_count;

    // Add all the contig.
    for (auto contig : genome.second->getMap()) {

      // Add all offsets.
      for (auto offset : contig.second->getMap()) {

        // if we find a homozygous variant then add it to the unphased.
        if (offset.second.size() == 1 and offset.second.front().second == vcf_ploidy) {

          std::shared_ptr<Variant> mutable_variant = std::const_pointer_cast<Variant>(offset.second.front().first);
          mutable_variant->updatePhaseId(ContigVariant::HAPLOID_HOMOLOGOUS_INDEX);   // Assign to the first (and only) homologous contig.
          genome_variant->addVariant(mutable_variant);
          ++homozygous_count;

        } else {  // the variant is heterozygous, we analyze the count statistics.

          ++heterozygous_count;

          size_t variant_index;
          if (analyseCountStatistics(offset.second, variant_index)) {

            std::shared_ptr<Variant> mutable_variant = std::const_pointer_cast<Variant>(offset.second.at(variant_index).first);
            mutable_variant->updatePhaseId(ContigVariant::HAPLOID_HOMOLOGOUS_INDEX);   // Assign to the first (and only) homologous contig.
            genome_variant->addVariant(mutable_variant);
            ++accepted_heterozygous_count;

          }

        }

      } // All offsets

    } // All contigs

  } // All genomes

  ExecEnv::log().info("Haploid Phasing; pre_phase variants: {}, genomes: {}, resultant population variants: {}, genomes: {}",
                      unphased_population_ptr->variantCount(), unphased_population_ptr->getMap().size(),
                      haploid_population->variantCount(), haploid_population->getMap().size());

  double percent_statistics = (static_cast<double>(statistics_count) / static_cast<double>(unphased_population_ptr->getMap().size())) * 100.0;
  double percent_mono = (static_cast<double>(mono_strain_count) / static_cast<double>(unphased_population_ptr->getMap().size())) * 100.0;

  ExecEnv::log().info("Haploid Phasing; Genomes with phasing statistics: {}%, Mono strain genomes: {}%", percent_statistics, percent_mono);

  double percent_homozygous = (static_cast<double>(homozygous_count) / static_cast<double>(homozygous_count + heterozygous_count)) * 100.0;
  double percent_accepted_heterozygous = (static_cast<double>(accepted_heterozygous_count) / static_cast<double>(heterozygous_count)) * 100.0;

  ExecEnv::log().info("Haploid Phasing; Homozygous variants: {}%, Accepted Heterozygous: {}% (proportion alternate >= {})",
                      percent_homozygous, percent_accepted_heterozygous, HETEROZYGOUS_PROPORTION_);

  return result;

}


bool kgl::GenomePhasing::analyseCountStatistics(const UnphasedVectorVariantCount& unphased_vector, size_t& phase_index) {

  phase_index = 0;
  for (auto variant : unphased_vector) {

    std::shared_ptr<const CountEvidence> read_evidence_ptr = std::dynamic_pointer_cast<const CountEvidence>(variant.first->evidence());

    if (not read_evidence_ptr) {

      ExecEnv::log().error("Haploid Phasing; Unexpected evidence ptr type for variant: {}",
                           variant.first->output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

    double proportion_alternate = static_cast<double>(read_evidence_ptr->altCount())
                                  / static_cast<double>(read_evidence_ptr->altCount() + read_evidence_ptr->refCount());

    if (proportion_alternate >= HETEROZYGOUS_PROPORTION_) {

      return true;

    }

    ++phase_index;

  }

  return false;

}
