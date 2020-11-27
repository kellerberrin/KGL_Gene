//
// Created by kellerberrin on 23/8/20.
//

#include <kgl_variant_factory_vcf_evidence_analysis.h>
#include "kgl_variant.h"
#include "kgl_filter.h"
#include "kgl_analysis_mutation_inbreed_locus.h"


namespace kgl = kellerberrin::genome;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//


std::vector<kgl::ContigOffset_t>  kgl::RetrieveLociiVector::getLociiVector(const LociiVectorArguments& arguments) {

  switch(arguments.lociiType()) {

    case LociiType::LOCII_FROM_COUNT:
      return getLociiCount(arguments.unphasedContigPtr(),
                           arguments.lowerOffset(),
                           arguments.lociiCount(),
                           arguments.lociiSpacing(),
                           arguments.superPopulation(),
                           arguments.minAlleleFrequency(),
                           arguments.maxAlleleFrequency(),
                           arguments.variantSource());

      case LociiType::LOCII_FROM_TO:
      return getLociiFromTo(arguments.unphasedContigPtr(),
                            arguments.lowerOffset(),
                            arguments.upperOffset(),
                            arguments.lociiSpacing(),
                            arguments.superPopulation(),
                            arguments.minAlleleFrequency(),
                            arguments.maxAlleleFrequency(),
                            arguments.variantSource());

  }

  return std::vector<kgl::ContigOffset_t>();

}



std::vector<kgl::AlleleFreqVector> kgl::RetrieveLociiVector::getAllelesFromTo( std::shared_ptr<const ContigVariant> unphased_contig_ptr,
                                                                             ContigOffset_t lower_offset,
                                                                             ContigOffset_t upper_offset,
                                                                             size_t locii_spacing,
                                                                             const std::string& super_population,
                                                                             double allele_frequency_min,
                                                                             double allele_frequency_max,
                                                                             VariantDatabaseSource variant_source) {

  std::vector<AlleleFreqVector> locii_vector;
  auto current_offset = unphased_contig_ptr->getMap().lower_bound(lower_offset);
  ContigOffset_t previous_offset{0};

  while (current_offset != unphased_contig_ptr->getMap().end()) {

    auto const&[offset, offset_ptr] = *current_offset;

    if (offset > upper_offset) {

      // No more locii.
      break;

    } else if ((offset >= previous_offset + locii_spacing) or previous_offset == 0) {

      OffsetVariantArray locus_variant_array = offset_ptr->getVariantArray();

      AlleleFreqVector allele_freq_vector(locus_variant_array, super_population, variant_source);

      if (not allele_freq_vector.checkValidAlleleVector()) {

        ++current_offset;  // next offset
        continue; // Next locus.

      }
      // Check for allele frequencies.
      double sum_frequencies = allele_freq_vector.minorAlleleFrequencies();
      size_t minor_allele_count = allele_freq_vector.alleleFrequencies().size();
      if (minor_allele_count == 0 or sum_frequencies == 0.0
          or sum_frequencies < allele_frequency_min or sum_frequencies > allele_frequency_max) {

        ++current_offset;  // next offset
        continue;

      }

      previous_offset = offset;
      locii_vector.push_back(allele_freq_vector);  // save offset for further use.

    } // if locii_spacing

    ++current_offset;  // next offset

  } // while

  return locii_vector;

}



std::vector<kgl::ContigOffset_t> kgl::RetrieveLociiVector::getLociiFromTo(std::shared_ptr<const ContigVariant> unphased_contig_ptr,
                                                                          ContigOffset_t lower_offset,
                                                                          ContigOffset_t upper_offset,
                                                                          size_t locii_spacing,
                                                                          const std::string& super_population,
                                                                          double allele_frequency_min,
                                                                          double allele_frequency_max,
                                                                          VariantDatabaseSource variant_source) {

  std::vector<kgl::ContigOffset_t> locii_vector;

  std::vector<AlleleFreqVector> freq_vector = getAllelesFromTo( unphased_contig_ptr,
                                                                lower_offset,
                                                                upper_offset,
                                                                locii_spacing,
                                                                super_population,
                                                                allele_frequency_min,
                                                                allele_frequency_max,
                                                                variant_source);

  for (auto const& alleles : freq_vector) {

    if (not alleles.alleleFrequencies().empty()) {

      locii_vector.push_back(alleles.alleleFrequencies().front().allele()->offset());

    } else {

      ExecEnv::log().error("RetrieveLociiVector::getLociiFromTo; empty allele vector");

    }

  }

  return locii_vector;

}




std::vector<kgl::AlleleFreqVector> kgl::RetrieveLociiVector::getAllelesCount( std::shared_ptr<const ContigVariant> unphased_contig_ptr,
                                                                              ContigOffset_t lower_offset,
                                                                              size_t locii_count,
                                                                              size_t locii_spacing,
                                                                              const std::string& super_population,
                                                                              double allele_frequency_min,
                                                                              double allele_frequency_max,
                                                                              VariantDatabaseSource variant_source) {

  std::vector<AlleleFreqVector> locii_vector;
  auto current_offset = unphased_contig_ptr->getMap().lower_bound(lower_offset);
  ContigOffset_t previous_offset{0};

  while (current_offset != unphased_contig_ptr->getMap().end()) {

    auto const&[offset, offset_ptr] = *current_offset;

    if (locii_vector.size() >= locii_count) {

      // No more locii.
      break;

    } else if ((offset >= previous_offset + locii_spacing) or previous_offset == 0) {

      OffsetVariantArray locus_variant_array = offset_ptr->getVariantArray();

      AlleleFreqVector allele_freq_vector(locus_variant_array, super_population, variant_source);

      if (not allele_freq_vector.checkValidAlleleVector()) {

        ++current_offset;  // next offset
        continue; // Next locus.

      }
      // Check for allele frequencies.
      double sum_frequencies = allele_freq_vector.minorAlleleFrequencies();
      size_t minor_allele_count = allele_freq_vector.alleleFrequencies().size();
      if (minor_allele_count == 0 or sum_frequencies == 0.0
          or sum_frequencies < allele_frequency_min or sum_frequencies > allele_frequency_max) {

        ++current_offset;  // next offset
        continue;

      }

      previous_offset = offset;
      locii_vector.push_back(allele_freq_vector);  // save offset for further use.

    } // if locii_spacing

    ++current_offset;  // next offset

  } // while

  return locii_vector;

}


std::vector<kgl::ContigOffset_t> kgl::RetrieveLociiVector::getLociiCount( std::shared_ptr<const ContigVariant> unphased_contig_ptr,
                                                                          ContigOffset_t lower_offset,
                                                                          size_t locii_count,
                                                                          size_t locii_spacing,
                                                                          const std::string& super_population,
                                                                          double allele_frequency_min,
                                                                          double allele_frequency_max,
                                                                          VariantDatabaseSource variant_source) {

  std::vector<kgl::ContigOffset_t> locii_vector;

  std::vector<AlleleFreqVector> freq_vector = getAllelesCount( unphased_contig_ptr,
                                                               lower_offset,
                                                               locii_count,
                                                               locii_spacing,
                                                               super_population,
                                                               allele_frequency_min,
                                                               allele_frequency_max,
                                                               variant_source);

  for (auto const& alleles : freq_vector) {

    if (not alleles.alleleFrequencies().empty()) {

      locii_vector.push_back(alleles.alleleFrequencies().front().allele()->offset());

    } else {

      ExecEnv::log().error("RetrieveLociiVector::getLociiCount; empty allele vector");

    }

  }

  return locii_vector;

}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//


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

  // Annotate the variant list with the super population frequency identifier
  std::shared_ptr<ContigVariant> locus_list(std::make_shared<ContigVariant>(super_population));
  size_t locii{0};

  if (unphased_ptr->getMap().size() == 1) {

    auto [genome_id, genome_ptr] = *(unphased_ptr->getMap().begin());

    auto contig_opt = genome_ptr->getContig(contig_id);

    if (contig_opt) {


      // Filter for SNP and 'PASS' variants.
      auto snp_contig_ptr = contig_opt.value()->filterVariants(AndFilter(SNPFilter(), PassFilter()));

      // Retrieve the locii that meet the conditions.
      LociiVectorArguments locii_args(snp_contig_ptr);
      locii_args.lociiType(LociiType::LOCII_FROM_TO);
      locii_args.lowerOffset(lower_offset);
      locii_args.upperOffset(upper_offset);
      locii_args.lociiSpacing(locii_spacing);
      locii_args.superPopulation(super_population);
      locii_args.minAlleleFrequency(min_af);
      locii_args.maxAlleleFrequency(max_af);

      std::vector<ContigOffset_t> locii_vector = RetrieveLociiVector::getLociiVector(locii_args);

      locii += locii_vector.size();

      for (auto locus : locii_vector) {

        auto variant_array_opt = snp_contig_ptr->findOffsetArray(locus);
        if (variant_array_opt) {

          for (auto const& variant : variant_array_opt.value()) {

            if (not locus_list->addVariant(variant)) {

              ExecEnv::log().error("InbreedingAnalysis::getLocusList, Could not add variant: {}",
                                   variant->output(',', VariantOutputIndex::START_0_BASED, false));

            }

          }

        } else {

          ExecEnv::log().error("InbreedingAnalysis::getLocusList, Could not add Contig: {}, No variant found at offset: {}",
                               snp_contig_ptr->contigId(), locus);

        }

      }

    } // if contig_opt

  } else {

    ExecEnv::log().error("InbreedingAnalysis::getLocusList, Unphased Population has {} Genomes, Expected 1",
                         unphased_ptr->getMap().size());

  }

  ExecEnv::log().info("Locus List for super population: {} contains: Locii: {}, SNPs: {}",
                      locus_list->contigId(), locii, locus_list->variantCount());

  return {super_population, locus_list};

}


