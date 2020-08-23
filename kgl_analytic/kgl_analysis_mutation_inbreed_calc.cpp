//
// Created by kellerberrin on 21/8/20.
//

#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"


#include <fstream>

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These functions calculate the inbreeding coefficient for an individual by looking at multiple locii.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the experimental Ritland multi-locus calculation developed in my document.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::InbreedingAnalysis::LocusResults
kgl::InbreedingAnalysis::multiLocus1( const GenomeId_t& genome_id,
                                    const std::shared_ptr<const DiploidContig>& contig_ptr,
                                    const std::string& super_population_field,
                                    const std::shared_ptr<const ContigVariant>& locus_list) const {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  LocusResults locus_results;
  locus_results.genome = genome_id;
  locus_results.inbred_allele_sum = 0.0;
  locus_results.homo_count = 0;
  locus_results.hetero_count = 0;
  locus_results.total_allele_count = 0;
  size_t loci_found = 0;

  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {
    // Join on the diploid contig.

    auto locus_variant_array = offset_ptr->getVariantArray();
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);

    if (diploid_variant_opt) {

      auto const& diploid_offset = diploid_variant_opt.value();

      bool found_flag = false;
      double allele_frequency = 0.0;
      for (auto const& locus_variant : locus_variant_array) {

        if (diploid_offset[0]->analogous(*locus_variant)) {
          // Found the matching locus allele.
          // Get the allele super population frequency
          auto [result, AF_value] = processFloatField(locus_variant, super_population_field);
          if (result and AF_value > 0.0 and AF_value < 1.0) {

            found_flag = true;
            ++loci_found;
            allele_frequency = AF_value;
            locus_results.total_allele_count += locus_variant_array.size();

          } // valid AF

          break; // No need to search further.

        } // Found locus variant

      } // for all locus variants

      if (found_flag) {

        double locus_inbreeding = 0.0;
        // Determine if the sample alternate allele is Hom/Het or Mixed.
        if (diploid_offset.size() == 1) {
          // The sample is alt allele heterozygous
          // Find the matching locus allele

          ++locus_results.hetero_count;

        } else if (diploid_offset.size() == 2) {

          if (diploid_offset[0]->homozygous(*diploid_offset[1])) {
            // The sample is alt allele homozygous
            // Find the matching locus allele

            // The implied possible allele count include the reference allele; thus (k - 1) = array.size()
            locus_inbreeding = ((1.0 / allele_frequency) - 1.0)/static_cast<double>(locus_variant_array.size());
            ++locus_results.homo_count;

          } else {
            // The sample has different alt alleles. Possible but unlikely.
            ExecEnv::log().info("InbreedingAnalysis::multiLocus1; Diploid genome: {} has two different non-ref alleles\n{}\n{}",
                                genome_id,
                                diploid_offset[0]->output(',',VariantOutputIndex::START_0_BASED, false),
                                diploid_offset[1]->output(',',VariantOutputIndex::START_0_BASED, false));

          }

        } else {

          // Sample has more than 2 alleles (should not happen if SNP filtered)
          ExecEnv::log().error("InbreedingAnalysis::multiLocus1; Diploid genome: {} has: {} SNPs at offset: {} contig: {}",
                               genome_id, diploid_offset.size(), offset, contig_ptr->contigId());
          continue;
        }

        locus_results.inbred_allele_sum += locus_inbreeding;

      }  // Sample allele found

    } // if sample locus

  } // for all loci.

  locus_results.inbred_allele_sum = locus_results.inbred_allele_sum / static_cast<double>(loci_found);

  ExecEnv::log().info("Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, Inbreeding: {}",
                      locus_results.genome, super_population_field, locus_results.hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum);

  return locus_results;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the published Ritland multi-locus calculation.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::InbreedingAnalysis::LocusResults
kgl::InbreedingAnalysis::processRitlandLocus(const GenomeId_t &genome_id,
                                           const std::shared_ptr<const DiploidContig>& contig_ptr,
                                           const std::string& super_population_field,
                                           const std::shared_ptr<const ContigVariant>& locus_list) const {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  LocusResults locus_results;
  locus_results.inbred_allele_sum = 0.0;
  locus_results.genome = genome_id;
  size_t sum_alternate_allele = 0;
  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {
    // Join on the diploid contig.

    auto locus_variant_array = offset_ptr->getVariantArray();
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);

    if (diploid_variant_opt) {

      sum_alternate_allele += locus_variant_array.size();
      ++locus_results.total_allele_count;
      // Determine if the sample alternate allele is Hom/Het or Mixed.
      auto const& diploid_offset = diploid_variant_opt.value();
      if (diploid_offset.size() == 1) {
        // The sample is alt allele heterozygous
        ++locus_results.hetero_count;
//        locus_results.inbred_allele_sum -= 1.0;

      } else if (diploid_offset.size() == 2) {

         if (diploid_offset[0]->homozygous(*diploid_offset[1])) {
          // The sample is alt allele homozygous
          // Find the matching locus allele
           ++locus_results.homo_count;

          for (auto const& locus_variant : locus_variant_array) {

            if (diploid_offset[0]->analogous(*locus_variant)) {
              // Found the matching locus allele.
              // Get the allele super population frequency
              auto [result, AF_value] = processFloatField(locus_variant, super_population_field);
              if (result and AF_value > 0.0 and AF_value < 1.0) {

                locus_results.inbred_allele_sum += (1.0 / AF_value);
//                locus_results.inbred_allele_sum += (1.0 / (1.0 - AF_value));
                locus_results.inbred_allele_sum -= 1.0;

              } // valid AF

              break; // No need to search further.

            } // Found locus variant

          } // for all locus variants

        } else {
          // The sample has different alt alleles.
          ExecEnv::log().info("InbreedingAnalysis::processRitlandLocus; Diploid genome: {} has two different non-ref alleles\n{}\n{}",
                              genome_id,
                              diploid_offset[0]->output(',',VariantOutputIndex::START_0_BASED, false),
                              diploid_offset[1]->output(',',VariantOutputIndex::START_0_BASED, false));

        }

      } else {

        ExecEnv::log().error("InbreedingAnalysis::processRitlandLocus; Diploid genome: {} has: {} SNPs at offset: {} contig: {}",
                             genome_id, diploid_offset.size(), offset, contig_ptr->contigId());
        continue;
      }

    }

  } // for all locus variants.

  locus_results.inbred_allele_sum = (sum_alternate_allele > 0 ? locus_results.inbred_allele_sum / static_cast<double>(sum_alternate_allele) : 0.0);

  ExecEnv::log().info("Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, Inbreeding: {}",
                      locus_results.genome, super_population_field, locus_results.hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum);

  return locus_results;

}

