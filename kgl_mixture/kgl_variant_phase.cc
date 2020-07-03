//
// Created by kellerberrin on 23/04/18.
//

#include "kgl_variant_phase.h"
#include "kgl_read_phasing.h"
#include "kgl_variant_evidence.h"
#include <fstream>


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Not Thread safe.
// This object accepts unphased variants and trivially phases them as haploid.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


// Just copy into a population object.
bool kgl::GenomePhasing::haploidPhasing(size_t vcf_ploidy,
                                        const std::shared_ptr<const UnphasedPopulation>& unphased_population_ptr,
                                        const std::shared_ptr<PhasedPopulation>& haploid_population)  {


  // Assumes a VCF ploidy of 2 and produces a haploid phased population.
  if (vcf_ploidy != 2) {

    ExecEnv::log().error("GenomePhasing::Haploid Phasing(); Expecting a VCF ploidy of 2, VCF ploidy is: {}", vcf_ploidy);
    return false;

  }

  size_t heterozygous_count = 0;
  size_t accepted_heterozygous_count = 0;
  size_t homozygous_count = 0;
  for (auto const& [genome_id, genome_ptr] : unphased_population_ptr->getMap()) {

    // Create the GenomeVariant object.
    std::optional<std::shared_ptr<GenomeVariant>> genome_variant_opt = haploid_population->getCreateGenome(genome_id);
    if (not genome_variant_opt) {

      ExecEnv::log().error("GenomePhasing::Haploid Phasing(); Unable to get/create genome: {} to haploid population", genome_id);
      return false;

    }

    // Add all the contig.
    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      // Add all offsets.
      for (auto const& [offset, variant_vector] : contig_ptr->getMap()) {

        // if we find a homozygous variant then add it to the unphased.
        if (variant_vector->getVariantArray().size() == 2 and variant_vector->getVariantArray()[0]->equivalent(*variant_vector->getVariantArray()[1])) {

          std::shared_ptr<Variant> mutable_variant = std::const_pointer_cast<Variant>(variant_vector->getVariantArray().front());
          mutable_variant->updatePhaseId(ContigVariant::HAPLOID_HOMOLOGOUS_INDEX);   // Assign to the first (and only) homologous contig.

          if (not genome_variant_opt.value()->addVariant(mutable_variant)) {

            ExecEnv::log().error("GenomePhasing::Haploid Phasing(); Haploid Genome: {}, Contig: {}, offset: {} Unable to add Homozygous variant",
                                  genome_id, contig_id, offset);

          }

          ++homozygous_count;

        } else if (variant_vector->getVariantArray().size() == 2 or variant_vector->getVariantArray().size() == 1) {  // the variant is heterozygous, we analyze the count statistics.

          ++heterozygous_count;

          size_t variant_index;
          if (analyseCountStatistics(*variant_vector, variant_index)) {

            std::shared_ptr<Variant> mutable_variant = std::const_pointer_cast<Variant>(variant_vector->getVariantArray()[variant_index]);
            mutable_variant->updatePhaseId(ContigVariant::HAPLOID_HOMOLOGOUS_INDEX);   // Assign to the first (and only) homologous contig.

            if (not genome_variant_opt.value()->addVariant(mutable_variant)) {

              ExecEnv::log().error("GenomePhasing::Haploid Phasing(); Haploid Genome: {}, Contig: {}, offset: {} Unable to add Heterozygous variant",
                                   genome_id, contig_id, offset);

            }

            ++accepted_heterozygous_count;

          }

        } else {

          ExecEnv::log().error("GenomePhasing::Haploid Phasing(); Haploid Genome: {}, Contig: {}, offset: {}, {} variants (> 2) at unique offset.",
                               genome_id, contig_id, offset, variant_vector->getVariantArray().size());
          for (auto const& variant_ptr : variant_vector->getVariantArray()) {

            ExecEnv::log().error("GenomePhasing::Haploid Phasing(); Variant: {}", variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false));

          }

        }

      } // All offsets

    } // All contigs

  } // All genomes

  ExecEnv::log().info("Haploid Phasing; pre_phase variants: {}, genomes: {}, resultant population variants: {}, genomes: {}",
                      unphased_population_ptr->variantCount(), unphased_population_ptr->getMap().size(),
                      haploid_population->variantCount(), haploid_population->getMap().size());

  if (heterozygous_count > 0) {

    double percent_homozygous = (static_cast<double>(homozygous_count) / static_cast<double>(homozygous_count + heterozygous_count)) * 100.0;
    double percent_accepted_heterozygous = (static_cast<double>(accepted_heterozygous_count) / static_cast<double>(heterozygous_count)) * 100.0;
    ExecEnv::log().info("Haploid Phasing; Homozygous variants: {}%, Accepted Heterozygous: {}% (proportion alternate >= {})",
                        percent_homozygous, percent_accepted_heterozygous, HETEROZYGOUS_PROPORTION_);

  } else {

    ExecEnv::log().info("Haploid Phasing; Homozygous variants: 100%");

  }


  return true;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Not Thread safe.
// This object accepts a file with phasing information (Pf3k only) and phases the variants for a haploid
// organism.
// Important - only samples with a Mixture coefficient of 1 are phased. Samples with
// multiple strains are ignored.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<kgl::UnphasedPopulation> kgl::GenomePhasing::filterClonal(const std::string& phase_file,
                                                                          std::shared_ptr<const UnphasedPopulation> unphased_population_ptr)  {

  PopulationId_t filtered_name = unphased_population_ptr->populationId() + std::string("_Clonal_Filter");
  std::shared_ptr<UnphasedPopulation> filtered_population(std::make_shared<UnphasedPopulation>(filtered_name));

  GenomeMixtureStatistics mixture_statistics;

  if (not mixture_statistics.readMixtureStatistics(phase_file)) {

    ExecEnv::log().info("GenomePhasing::filterClonal; cannot read phasing file: {}", phase_file);
    return filtered_population;

  }

  for (auto genome : unphased_population_ptr->getMap()) {


    // Create the GenomeVariant object.

    MixtureStatistics genome_statistics;
    if (not mixture_statistics.getMixtureStatistics(genome.first, genome_statistics)) {

      ExecEnv::log().info("GenomePhasing::filterClonal; No phasing statistics available for genome: {}", genome.first);
      continue;

    }

    if (genome_statistics.first > 1) {

      continue;

    }

    // Is a mono (clonal) strain so create the phased genome object.

    if (not filtered_population->addGenome(genome.second->deepCopy())) {

      ExecEnv::log().error("GenomePhasing::filterClonal; Unable to add clonal genome: {} to population", genome.first);

    }

  }

  ExecEnv::log().info("GenomePhasing::filterClonal; Total Genomes: {}, Clonal Genomes: {}",
                      unphased_population_ptr->getMap().size(), filtered_population->getMap().size());

  return filtered_population;

}


bool kgl::GenomePhasing::analyseCountStatistics(const VirtualContigOffset& unphased_vector, size_t& phase_index) {

  phase_index = 0;
  double proportion_alternate = 0.0;
  for (size_t i = 0; i < unphased_vector.getVariantArray().size(); ++i) {

    std::optional<std::shared_ptr<const FormatData>> read_evidence_opt = unphased_vector.getVariantArray()[i]->evidence().formatData();

    if (not read_evidence_opt) {

      ExecEnv::log().error("Haploid Phasing; Unexpected evidence ptr type for variant: {}",
                           unphased_vector.getVariantArray()[i]->output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

    auto total_counts = read_evidence_opt.value()->DPCount();
    if (total_counts == 0) continue;

    double proportion_variant = static_cast<double>(read_evidence_opt.value()->altCount()) / static_cast<double>(total_counts);

    if (proportion_variant > proportion_alternate) {

      proportion_alternate = proportion_variant;
      phase_index = i;

    }

  }

  return proportion_alternate >= HETEROZYGOUS_PROPORTION_;

}
