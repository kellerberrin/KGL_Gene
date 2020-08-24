//
// Created by kellerberrin on 23/8/20.
//

//
// Created by kellerberrin on 23/8/20.
//

#include <kgl_variant_factory_vcf_evidence_analysis.h>
#include "kgl_variant.h"
#include "kgl_filter.h"
#include "kgl_analysis_mutation_inbreed_aux.h"


namespace kgl = kellerberrin::genome;


std::pair<bool, double> kgl::InbreedSampling::processFloatField( const Variant& variant,
                                                                 const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(variant, field_name);

  if (field_opt) {

    std::vector<double> field_vec = InfoEvidenceAnalysis::varianttoFloats(field_opt.value());

    if (field_vec.size() == 1) {

      return {true, field_vec.front() };

    } else if (field_vec.size() == 0) {

      // Missing value
      return {false, 0.0};

    } else {

      std::string vector_str;
      for (auto const& str : field_vec) {

        vector_str += str;
        vector_str += ";";

      }

      ExecEnv::log().error("InbreedingAnalysis::processField, Field: {} expected vector size 1, get vector size: {}, vector: {}",
                           field_name, field_vec.size(), vector_str);
      return {false, 0.0};

    }

  } else {

    ExecEnv::log().error("InbreedingAnalysis::processField, Field: {} not found for Variant: {}",
                         field_name, variant.output(',',VariantOutputIndex::START_0_BASED, false));

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


kgl::LocusMap kgl::InbreedSampling::getPopulationLocus(const UnphasedPopulation& unphased_population,
                                                       const ContigId_t& contig_id) {

  LocusMap locus_map;

  // Get an allele sampling locus for each of the defined superpopulations.
  locus_map[SUPER_POP_AFR_GNOMAD_.first] = getLocusList(unphased_population,
                                                        contig_id,
                                                        VARIANT_SPACING_,
                                                        SUPER_POP_AFR_GNOMAD_.second,
                                                        MIN_MAF,
                                                        MAX_MAF);

  locus_map[SUPER_POP_AMR_GNOMAD_.first] = getLocusList(unphased_population,
                                                        contig_id,
                                                        VARIANT_SPACING_,
                                                        SUPER_POP_AMR_GNOMAD_.second,
                                                        MIN_MAF,
                                                        MAX_MAF);

  locus_map[SUPER_POP_EAS_GNOMAD_.first] = getLocusList(unphased_population,
                                                        contig_id,
                                                        VARIANT_SPACING_,
                                                        SUPER_POP_EAS_GNOMAD_.second,
                                                        MIN_MAF,
                                                        MAX_MAF);

  locus_map[SUPER_POP_EUR_GNOMAD_.first] = getLocusList(unphased_population,
                                                        contig_id,
                                                        VARIANT_SPACING_,
                                                        SUPER_POP_EUR_GNOMAD_.second,
                                                        MIN_MAF,
                                                        MAX_MAF);

  locus_map[SUPER_POP_SAS_GNOMAD_.first] = getLocusList(unphased_population,
                                                        contig_id,
                                                        VARIANT_SPACING_,
                                                        SUPER_POP_SAS_GNOMAD_.second,
                                                        MIN_MAF,
                                                        MAX_MAF);

  return locus_map;

}



// Get a list of hom/het SNPs with a spcified spacing to minimise linkage dis-equilibrium
// and at a specified frequency for the super population. Used as a template for calculating
// the inbreeding coefficient and sample relatedness
std::shared_ptr<const kgl::ContigVariant> kgl::InbreedSampling::getLocusList( const UnphasedPopulation& unphased_population,
                                                                              const ContigId_t& contig_id,
                                                                              ContigOffset_t spacing,
                                                                              const std::string& super_pop_freq,
                                                                              double min_frequency,
                                                                              double max_frequency) {

  // Annotate the variant list with the super population frequency identifier
  std::shared_ptr<ContigVariant> locus_list(std::make_shared<ContigVariant>(super_pop_freq));

  if (unphased_population.getMap().size() == 1) {

    auto [genome_id, genome_ptr] = *(unphased_population.getMap().begin());

    auto contig_opt = genome_ptr->getContig(contig_id);

    if (contig_opt) {

      // Filter for SNP.
      auto snp_contig_ptr = contig_opt.value()->filterVariants(SNPFilter());
      // Filter for maximum and minimum AF frequency
      snp_contig_ptr = snp_contig_ptr->filterVariants(AndFilter(InfoGEQFloatFilter(super_pop_freq, min_frequency),
                                                                NotFilter(InfoGEQFloatFilter(super_pop_freq, max_frequency))));

      ContigOffset_t previous_offset{0};
      for (auto const& [offset, offset_ptr] : snp_contig_ptr->getMap()) {

        if (offset >= previous_offset + spacing) {

          OffsetVariantArray variant_array = offset_ptr->getVariantArray();

          for (auto const& variant_ptr : variant_array) {

            if (not locus_list->addVariant(variant_ptr)) {

              ExecEnv::log().error("InbreedingAnalysis::getLocusList, Could not add variant: {}",
                                   variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

            }

          } // for variant array

        } // if spacing

      } // for offset

    } // if contig_opt

  } else {

    ExecEnv::log().error("InbreedingAnalysis::getLocusList, Unphased Population has {} Genomes, Expected 1",
                         unphased_population.getMap().size());

  }

  if (locus_list->variantCount() > 0) {

    ExecEnv::log().info("Locus List for super population: {} contains: {} SNPs", locus_list->contigId(), locus_list->variantCount());

  }

  return locus_list;

}


