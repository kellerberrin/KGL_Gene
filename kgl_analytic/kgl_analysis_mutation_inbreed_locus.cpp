//
// Created by kellerberrin on 23/8/20.
//

#include <kgl_variant_factory_vcf_evidence_analysis.h>
#include "kgl_variant.h"
#include "kgl_filter.h"
#include "kgl_analysis_mutation_inbreed_locus.h"
#include "kgl_analysis_mutation_inbreed_calc.h"


namespace kgl = kellerberrin::genome;



kgl::ContigLocusMap kgl::InbreedSampling::getPopulationLocusMap(  std::shared_ptr<const UnphasedPopulation> population_ptr,
                                                                  double min_af,
                                                                  double max_af,
                                                                  ContigOffset_t locii_spacing,
                                                                  ContigOffset_t upper_offset,
                                                                  ContigOffset_t lower_offset) {

  ContigLocusMap contig_locus_map;
  // We assume that the unphased population has only 1 genome.
  if (population_ptr->getMap().size() == 1) {

    // Retrieve the genome ptr.
    auto& [genome_id, genome_ptr] = *(population_ptr->getMap().begin());
    // Generate contig locii maps
    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      contig_locus_map[contig_id] = getPopulationLocus( population_ptr,
                                                        contig_id,
                                                        min_af,
                                                        max_af,
                                                        locii_spacing,
                                                        upper_offset,
                                                        lower_offset);


      ExecEnv::log().info( "Generated locus maps for contig: {}", contig_id);

    }

    return contig_locus_map;

  } else {
  // If not 1 genome, issue an error message and return an empty map;

    ExecEnv::log().error( "InbreedSampling::getPopulationLocusMap, The Gnomad population has: {} genomes, expected 1.",
                          population_ptr->getMap().size());

    return contig_locus_map;

  }

}


kgl::LocusMap kgl::InbreedSampling::getPopulationLocus(std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                       const ContigId_t& contig_id,
                                                       double min_af,
                                                       double max_af,
                                                       ContigOffset_t locii_spacing,
                                                       ContigOffset_t upper_offset,
                                                       ContigOffset_t lower_offset) {


  ThreadPool threadpool(VariantDatabaseRead::superPopulations().size());
  std::vector<std::future<LocusReturnPair>> futures_vec;

  for (auto const& super_pop : VariantDatabaseRead::superPopulations()) {

    auto return_future = threadpool.enqueueTask(&InbreedSampling::getLocusList,
                                                unphased_ptr,
                                                contig_id,
                                                super_pop,
                                                min_af,
                                                max_af,
                                                locii_spacing,
                                                upper_offset,
                                                lower_offset);

    futures_vec.push_back(std::move(return_future));

  }

  LocusMap locus_map;
  for (auto& future : futures_vec) {

    auto [ super_population, locus_list] = future.get();

    locus_map[super_population] = locus_list;

  }

  return locus_map;

}



// Get a list of hom/het SNPs with a specified locii_spacing to minimise linkage dis-equilibrium
// and at a specified frequency for the super population. Used as a template for calculating
// the inbreeding coefficient and sample relatedness
kgl::InbreedSampling::LocusReturnPair kgl::InbreedSampling::getLocusList( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                                          const ContigId_t& contig_id,
                                                                          const std::string& super_population,
                                                                          double min_af,
                                                                          double max_af,
                                                                          ContigOffset_t locii_spacing,
                                                                          ContigOffset_t upper_offset,
                                                                          ContigOffset_t lower_offset) {

  size_t locii{0};
  size_t variant_duplicates{0};
  static std::mutex log_mutex;  // Multi-thread logging.
  // Annotate the variant list with the super population frequency identifier
  std::shared_ptr<ContigVariant> locus_list(std::make_shared<ContigVariant>(super_population));

  if (unphased_ptr->getMap().size() == 1) {

    auto [genome_id, genome_ptr] = *(unphased_ptr->getMap().begin());

    auto contig_opt = genome_ptr->getContig(contig_id);

    if (contig_opt) {

      // Filter for SNP and 'PASS' variants.
      auto snp_contig_ptr = contig_opt.value()->filterVariants(AndFilter(SNPFilter(), PassFilter()));

      ContigOffset_t previous_offset = lower_offset;
      for (auto const& [offset, offset_ptr] : snp_contig_ptr->getMap()) {

        if (offset > upper_offset) {

          // No more locii.
          break;

        } else if (offset >= previous_offset + locii_spacing) {

          OffsetVariantArray locus_variant_array = offset_ptr->getVariantArray();

          AlleleFreqVector allele_freq_vector(locus_variant_array, super_population);


          if (not allele_freq_vector.checkValidAlleleVector()) {

            continue; // Next locus.

          }
          // Check for allele frequencies.
          double sum_frequencies = allele_freq_vector.minorAlleleFrequencies();
          size_t minor_allele_count = allele_freq_vector.alleleFrequencies().size();
          if (minor_allele_count == 0 or sum_frequencies == 0.0
              or sum_frequencies < min_af or sum_frequencies > max_af) {

            continue;

          }

          for (auto const& allele : allele_freq_vector.alleleFrequencies()) {

            if (not locus_list->addVariant(allele.allele())) {

              ExecEnv::log().error("InbreedingAnalysis::getLocusList, Could not add variant: {}",
                                   allele.allele()->output(',', VariantOutputIndex::START_0_BASED, false));

            }

          } // for variant array

          ++locii;
          previous_offset = offset;

        } // if locii_spacing

      } // for offset

    } // if contig_opt

  } else {

    ExecEnv::log().error("InbreedingAnalysis::getLocusList, Unphased Population has {} Genomes, Expected 1",
                         unphased_ptr->getMap().size());

  }

  if (locus_list->variantCount() > 0) {

    ExecEnv::log().info("Locus List for super population: {} contains: Locii: {}, SNPs: {}, rejected duplicates: {}",
                        locus_list->contigId(), locii, locus_list->variantCount(), variant_duplicates);

  }

  return {super_population, locus_list};

}


