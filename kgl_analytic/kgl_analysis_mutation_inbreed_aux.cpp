//
// Created by kellerberrin on 23/8/20.
//

#include <kgl_variant_factory_vcf_evidence_analysis.h>
#include "kgl_variant.h"
#include "kgl_filter.h"
#include "kgl_analysis_mutation_inbreed_aux.h"
#include "kgl_analysis_mutation_inbreed_calc.h"


namespace kgl = kellerberrin::genome;


std::pair<bool, double> kgl::InbreedSampling::processFloatField( const Variant& variant,
                                                                 const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(variant, field_name);

  if (field_opt) {

    std::vector<double> field_vec = InfoEvidenceAnalysis::varianttoFloats(field_opt.value());

    if (field_vec.size() == 1) {

      return {true, field_vec.front() };

    } else if (field_vec.size() == 0) {

      // Missing value, this is OK, means the field exists but not defined for this variant.
      return {false, 0.0};

    } else {

      std::string vector_str;
      for (auto const& str : field_vec) {

        vector_str += str;
        vector_str += ";";

      }

      ExecEnv::log().warn("InbreedingAnalysis::processFloatField, Field: {} expected vector size 1, get vector size: {}, vector: {}, Variant: {}",
                           field_name, field_vec.size(), vector_str, variant.output(',', VariantOutputIndex::START_0_BASED, false));
      return {false, 0.0};

    }

  } else {

    ExecEnv::log().warn("InbreedingAnalysis::processFloatField, Field: {} Not found for variant: {}",
                         field_name, variant.output(',', VariantOutputIndex::START_0_BASED, false));
    return {false, 0.0};

  }

}


std::string kgl::InbreedSampling::lookupSuperPopulationField(const std::string& super_population) {

  if (super_population == SUPER_POP_AFR_GNOMAD_.first) {

    return SUPER_POP_AFR_GNOMAD_.second;

  } else if (super_population == SUPER_POP_AMR_GNOMAD_.first) {

    return SUPER_POP_AMR_GNOMAD_.second;

  } else if (super_population == SUPER_POP_EAS_GNOMAD_.first) {

    return SUPER_POP_EAS_GNOMAD_.second;

  } else if (super_population == SUPER_POP_EUR_GNOMAD_.first) {

    return SUPER_POP_EUR_GNOMAD_.second;

  } else if (super_population == SUPER_POP_SAS_GNOMAD_.first) {

    return SUPER_POP_SAS_GNOMAD_.second;

  } else  {

    ExecEnv::log().error("MutationAnalysis::lookupSuperPopulationField; Unknown Super Population: {}", super_population);
    return SUPER_POP_SAS_GNOMAD_.second;

  }

}


std::string kgl::InbreedSampling::inverseSuperPopulationField(const std::string& super_population) {

  if (super_population == SUPER_POP_AFR_GNOMAD_.second) {

    return SUPER_POP_AFR_GNOMAD_.first + std::string("_AF");

  } else if (super_population == SUPER_POP_AMR_GNOMAD_.second) {

    return SUPER_POP_AMR_GNOMAD_.first + std::string("_AF");

  } else if (super_population == SUPER_POP_EAS_GNOMAD_.second) {

    return SUPER_POP_EAS_GNOMAD_.first + std::string("_AF");

  } else if (super_population == SUPER_POP_EUR_GNOMAD_.second) {

    return SUPER_POP_EUR_GNOMAD_.first + std::string("_AF");

  } else if (super_population == SUPER_POP_SAS_GNOMAD_.second) {

    return SUPER_POP_SAS_GNOMAD_.first + std::string("_AF");

  } else  {

    ExecEnv::log().error("MutationAnalysis::lookupSuperPopulationField; Unknown Super Population: {}", super_population);
    return SUPER_POP_SAS_GNOMAD_.first + std::string("_AF");

  }

}


kgl::ContigLocusMap kgl::InbreedSampling::getPopulationLocusMap(  std::shared_ptr<const UnphasedPopulation> population_ptr,
                                                                  double min_af,
                                                                  double max_af,
                                                                  ContigOffset_t locii_spacing) {

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
                                                        locii_spacing);


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
                                                       ContigOffset_t locii_spacing) {


  ThreadPool threadpool(super_pop_vec.size());
  std::vector<std::future<LocusReturnPair>> futures_vec;

  for (auto const& super_pop : super_pop_vec) {

    auto return_future = threadpool.enqueueTask(&InbreedSampling::getLocusList,
                                                unphased_ptr,
                                                contig_id,
                                                locii_spacing,
                                                super_pop,
                                                min_af,
                                                max_af);

    futures_vec.push_back(std::move(return_future));

  }

  LocusMap locus_map;
  for (auto& future : futures_vec) {

    auto [ super_population, locus_list] = future.get();

    locus_map[super_population] = locus_list;

  }

  return locus_map;

}



// Get a list of hom/het SNPs with a specified spacing to minimise linkage dis-equilibrium
// and at a specified frequency for the super population. Used as a template for calculating
// the inbreeding coefficient and sample relatedness
kgl::InbreedSampling::LocusReturnPair kgl::InbreedSampling::getLocusList( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                                          const ContigId_t& contig_id,
                                                                          ContigOffset_t spacing,
                                                                          const SuperPopPair& super_pop_pair,
                                                                          double min_frequency,
                                                                          double max_frequency) {

  size_t locii{0};
  size_t variant_duplicates{0};
  static std::mutex log_mutex;  // Multi-thread logging.
  // Annotate the variant list with the super population frequency identifier
  std::shared_ptr<ContigVariant> locus_list(std::make_shared<ContigVariant>(super_pop_pair.second));

  if (unphased_ptr->getMap().size() == 1) {

    auto [genome_id, genome_ptr] = *(unphased_ptr->getMap().begin());

    auto contig_opt = genome_ptr->getContig(contig_id);

    if (contig_opt) {

      // Filter for SNP and 'PASS' variants.
      auto snp_contig_ptr = contig_opt.value()->filterVariants(AndFilter(SNPFilter(), PassFilter()));

      ContigOffset_t previous_offset{0};
      for (auto const& [offset, offset_ptr] : snp_contig_ptr->getMap()) {

        if (offset >= previous_offset + spacing) {

          OffsetVariantArray locus_variant_array = offset_ptr->getVariantArray();

          AlleleFreqVector allele_freq_vector(locus_variant_array, super_pop_pair.second);


          if (not allele_freq_vector.checkValidAlleleVector()) {

            continue; // Next locus.

          }
          // Check for allele frequencies.
          double sum_frequencies = allele_freq_vector.minorAlleleFrequencies();
          size_t minor_allele_count = allele_freq_vector.alleleFrequencies().size();
          if (minor_allele_count == 0 or sum_frequencies == 0.0
              or sum_frequencies < min_frequency or sum_frequencies > max_frequency) {

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

        } // if spacing

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

  return {super_pop_pair.first, locus_list};

}


