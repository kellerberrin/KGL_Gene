//
// Created by kellerberrin on 13/8/20.
//


#include "kel_thread_pool.h"
#include "kgl_analysis_mutation.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"


#include <fstream>

namespace kgl = kellerberrin::genome;


// Calculate the HetHom Ratio
bool kgl::MutationAnalysis::hetHomRatioLocus(const std::shared_ptr<const DiploidPopulation>& population) const {

  for (auto const& [genome_contig_id, genome_contig_ptr] : genome_GRCh38_->getMap()) {

    // Generate the super population allele locus lists
    std::map<std::string, std::shared_ptr<const ContigVariant>> locus_map;
    locus_map[SUPER_POP_AFR_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_AFR_GNOMAD_.second, MIN_MAF, MAX_MAF);
    locus_map[SUPER_POP_AMR_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_AMR_GNOMAD_.second, MIN_MAF, MAX_MAF);
    locus_map[SUPER_POP_EAS_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_EAS_GNOMAD_.second, MIN_MAF, MAX_MAF);
    locus_map[SUPER_POP_EUR_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_EUR_GNOMAD_.second, MIN_MAF, MAX_MAF);
    locus_map[SUPER_POP_SAS_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_SAS_GNOMAD_.second, MIN_MAF, MAX_MAF);

    // Use a thread pool to calculate inbreeding and relatedness.
    ThreadPool thread_pool;
    std::vector<std::future<LocusResults>> future_vector;
    for (auto const&[genome_id, genome_ptr] : population->getMap()) {

      auto contig_opt = genome_ptr->getContig(genome_contig_id);

      if (contig_opt) {

        auto result = ped_data_->getMap().find(genome_id);
        if (result == ped_data_->getMap().end()) {

          ExecEnv::log().error("MutationAnalysis::hetHomRatioLocu, Genome sample: {} does not have a PED record", genome_id);
          continue;

        }
        auto const& [sample_id, ped_record] = *result;
        auto locus_result = locus_map.find(ped_record.superPopulation());
        if (locus_result == locus_map.end()) {

          ExecEnv::log().error("MutationAnalysis::hetHomRatioLocu, Locus set not found for super population: {}", ped_record.superPopulation());
          continue;

        }
        auto const& [super_pop_id, locus_list] = *locus_result;

        std::future<LocusResults> future = thread_pool.enqueueTask(&MutationAnalysis::multiLocus1,
                                                                   this,
                                                                   genome_id,
                                                                   contig_opt.value(),
                                                                   lookupSuperPopulationField(super_pop_id),
                                                                   locus_list);
        future_vector.push_back(std::move(future));

      }

    }

    // Retrieve the thread results into a map.
    std::map<GenomeId_t, LocusResults> genome_results_map;
    for (auto &future : future_vector) {

      auto locus_results = future.get();
      ExecEnv::log().info("Processed Diploid genome: {} for inbreeding and relatedness", locus_results.genome);
      genome_results_map[locus_results.genome] = locus_results;


    }

    if (not genome_results_map.empty()) {

      writeResults(genome_contig_id, genome_results_map);

    }

  }

  return true;

}


// Calculate the HetHom Ratio
bool kgl::MutationAnalysis::syntheticInbreeding() const {

  for (auto const& [genome_contig_id, genome_contig_ptr] : genome_GRCh38_->getMap()) {

    // Generate the super population allele locus lists
    std::map<std::string, std::shared_ptr<const ContigVariant>> locus_map;
    locus_map[SUPER_POP_AFR_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_AFR_GNOMAD_.second, MIN_MAF, MAX_MAF);
    locus_map[SUPER_POP_AMR_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_AMR_GNOMAD_.second, MIN_MAF, MAX_MAF);
    locus_map[SUPER_POP_EAS_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_EAS_GNOMAD_.second, MIN_MAF, MAX_MAF);
    locus_map[SUPER_POP_EUR_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_EUR_GNOMAD_.second, MIN_MAF, MAX_MAF);
    locus_map[SUPER_POP_SAS_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_SAS_GNOMAD_.second, MIN_MAF, MAX_MAF);

    // Use a thread pool to calculate inbreeding and relatedness.
    for (auto const&[super_pop_id, locus_list] : locus_map) {

      std::shared_ptr<const DiploidPopulation> population = generateSyntheticPopulation(MIN_INBREEDING_COEFICIENT,
                                                                                        MAX_INBREEDING_COEFICIENT,
                                                                                        STEP_INBREEDING_COEFICIENT,
                                                                                        super_pop_id,
                                                                                        locus_list);

      ThreadPool thread_pool;
      std::vector<std::future<LocusResults>> future_vector;
      for (auto const&[genome_id, genome_ptr] : population->getMap()) {

        auto contig_opt = genome_ptr->getContig(genome_contig_id);

        if (contig_opt) {


          std::future<LocusResults> future = thread_pool.enqueueTask(&MutationAnalysis::multiLocus1,
                                                                     this,
                                                                     genome_id,
                                                                     contig_opt.value(),
                                                                     lookupSuperPopulationField(super_pop_id),
                                                                     locus_list);
          future_vector.push_back(std::move(future));

        }

      }

      // Retrieve the thread results into a map.
      std::map<GenomeId_t, LocusResults> genome_results_map;
      for (auto &future : future_vector) {

        auto locus_results = future.get();
        ExecEnv::log().info("Processed Diploid genome: {} for inbreeding and relatedness", locus_results.genome);
        genome_results_map[locus_results.genome] = locus_results;


      }

      if (not genome_results_map.empty()) {

        syntheticResults(genome_contig_id, genome_results_map);

      }

    } // for locus

  } // for contig.

  return true;

}



bool kgl::MutationAnalysis::syntheticResults( const ContigId_t& contig_id,
                                              const std::map<GenomeId_t, LocusResults>& genome_results_map) const {

  // Append the results.
  std::ofstream outfile;
  outfile.open(output_file_name_, std::ofstream::out |  std::ofstream::app);

  outfile << contig_id << DELIMITER_ << "HetCount" << DELIMITER_ << "HomCount"
          << DELIMITER_ << "Het/Hom"<< DELIMITER_ << "TotalLoci" << DELIMITER_ << "Ritland\n";

  for (auto const& [genome_id, locus_results] : genome_results_map) {


    double het_hom_ratio = (locus_results.homo_count > 0 ? static_cast<double>(locus_results.hetero_count)/static_cast<double>(locus_results.homo_count) : 0.0);

    outfile << genome_id << DELIMITER_;
    outfile << locus_results.hetero_count << DELIMITER_;
    outfile << locus_results.homo_count << DELIMITER_;
    outfile << het_hom_ratio << DELIMITER_;
    outfile << locus_results.total_allele_count << DELIMITER_;
    outfile << locus_results.inbred_allele_sum << DELIMITER_;
    outfile << '\n';

  }

  outfile.flush();

  return true;

}


bool kgl::MutationAnalysis::writeResults( const ContigId_t& contig_id,
                                          const std::map<GenomeId_t, LocusResults>& genome_results_map) const {

  // Append the results.
  std::ofstream outfile;
  outfile.open(output_file_name_, std::ofstream::out |  std::ofstream::app);

  outfile << contig_id << DELIMITER_ << "Population" << DELIMITER_
          << "Description" << DELIMITER_ << "SuperPopulation" << DELIMITER_ << "Description" << DELIMITER_
          << "HetCount" << DELIMITER_ << "HomCount" << DELIMITER_ << "Het/Hom"<< DELIMITER_
          << "TotalLoci" << DELIMITER_ << "Ritland\n";

  for (auto const& [genome_id, locus_results] : genome_results_map) {

    auto result = ped_data_->getMap().find(genome_id);

    if (result == ped_data_->getMap().end()) {

      ExecEnv::log().error("MutationAnalysis::hetHomRatio, Genome sample: {} does not have a PED record", genome_id);
      continue;

    }

    auto const& [sample_id, ped_record] = *result;

    double het_hom_ratio = (locus_results.homo_count > 0 ? static_cast<double>(locus_results.hetero_count)/static_cast<double>(locus_results.homo_count) : 0.0);

    outfile << genome_id << DELIMITER_;
    outfile << ped_record.population() << DELIMITER_;
    outfile << ped_record.populationDescription() << DELIMITER_;
    outfile << ped_record.superPopulation() << DELIMITER_;
    outfile << ped_record.superDescription() << DELIMITER_;
    outfile << locus_results.hetero_count << DELIMITER_;
    outfile << locus_results.homo_count << DELIMITER_;
    outfile << het_hom_ratio << DELIMITER_;
    outfile << locus_results.total_allele_count << DELIMITER_;
    outfile << locus_results.inbred_allele_sum << DELIMITER_;
    outfile << '\n';

  }

  outfile.flush();

  return true;

}




// Experimental Ritland multi-locus
kgl::MutationAnalysis::LocusResults
kgl::MutationAnalysis::multiLocus1( const GenomeId_t& genome_id,
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
            ExecEnv::log().info("MutationAnalysis::processLocusContig; Diploid genome: {} has two different non-ref alleles\n{}\n{}",
                                genome_id,
                                diploid_offset[0]->output(',',VariantOutputIndex::START_0_BASED, false),
                                diploid_offset[1]->output(',',VariantOutputIndex::START_0_BASED, false));

          }

        } else {

        // Sample has more than 2 alleles (should not happen if SNP filtered)
          ExecEnv::log().error("MutationAnalysis::processLocusContig; Diploid genome: {} has: {} SNPs at offset: {} contig: {}",
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


kgl::MutationAnalysis::LocusResults
kgl::MutationAnalysis::processLocusContig( const GenomeId_t& genome_id,
                                           const std::shared_ptr<const DiploidContig>& contig_ptr,
                                           const std::string& super_population_field,
                                           const std::shared_ptr<const ContigVariant>& locus_list) const {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  LocusResults locus_results;
  locus_results.inbred_allele_sum = 0.0;
  locus_results.genome = genome_id;
  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {
    // Join on the diploid contig.

    auto locus_variant_array = offset_ptr->getVariantArray();
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);

    if (diploid_variant_opt) {
      // Determine if the sample alternate allele is Hom/Het or Mixed.
      auto const& diploid_offset = diploid_variant_opt.value();
      if (diploid_offset.size() == 1) {
        // The sample is alt allele heterozygous
        // Find the matching locus allele
        for (auto const& locus_variant : locus_variant_array) {

          if (diploid_offset[0]->analogous(*locus_variant)) {
            // Found the matching locus allele.
            // Get the allele super population frequency
            auto [result, AF_value] = processFloatField(locus_variant, super_population_field);
            if (result and AF_value > 0.0 and AF_value < 1.0) {

              ++locus_results.total_allele_count;
              locus_results.inbred_allele_sum += -1.0 * AF_value;
              locus_results.inbred_allele_sum += (1.0 / (1.0 - AF_value)) - (1.0 - AF_value);

            } // Valid AF

            ++locus_results.hetero_count;
            break; // No need to search further.

          } // Found locus variant

        } // for all locus variants.

      } else if (diploid_offset.size() == 2) {

        if (diploid_offset[0]->analogous(*diploid_offset[1])) {
          // The sample is alt allele homozygous
          // Find the matching locus allele
          for (auto const& locus_variant : locus_variant_array) {

            if (diploid_offset[0]->analogous(*locus_variant)) {
              // Found the matching locus allele.
              // Get the allele super population frequency
              auto [result, AF_value] = processFloatField(locus_variant, super_population_field);
              if (result and AF_value > 0.0 and AF_value < 1.0 and false) {

                ++locus_results.total_allele_count;
                locus_results.inbred_allele_sum += (1.0 / AF_value) - AF_value;
                locus_results.inbred_allele_sum += (1.0 / (1.0 - AF_value)) - (1.0 - AF_value);

              } // valid AF

              ++locus_results.homo_count;
              break; // No need to search further.

            } // Found locus variant

          } // for all locus variants

        } else {
          // The sample has different alt alleles.
          ExecEnv::log().info("MutationAnalysis::processLocusContig; Diploid genome: {} has two different non-ref alleles\n{}\n{}",
                              genome_id,
                              diploid_offset[0]->output(',',VariantOutputIndex::START_0_BASED, false),
                              diploid_offset[1]->output(',',VariantOutputIndex::START_0_BASED, false));

        }

      } else {

        ExecEnv::log().error("MutationAnalysis::processLocusContig; Diploid genome: {} has: {} SNPs at offset: {} contig: {}",
                             genome_id, diploid_offset.size(), offset, contig_ptr->contigId());
        continue;
      }

    } else {
      // The sample has a homozygous reference allele at this location.
      double sum_AF_value{0.0};
      for (auto const& locus_variant : locus_variant_array) {

          // Found the matching locus allele.
          // Get the minor allele super population frequency
        auto [result, AF_value] = processFloatField(locus_variant, super_population_field);
        if (result and AF_value > 0.0) {

          ++locus_results.total_allele_count;
          locus_results.inbred_allele_sum += (1.0 / AF_value) - AF_value;
          sum_AF_value += AF_value;

        } // valid AF

      } // for all locus variants

      if (sum_AF_value > 0.0 and sum_AF_value < 1.0) {

        locus_results.inbred_allele_sum += (1.0 / (1.0 - sum_AF_value)) - (1.0 - sum_AF_value);

      }

    }

  } // for all locus variants.

  locus_results.inbred_allele_sum = (locus_results.total_allele_count > 0 ? locus_results.inbred_allele_sum / static_cast<double>(locus_results.total_allele_count) : 0.0);

  ExecEnv::log().info("Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, Inbreeding: {}",
                      locus_results.genome, super_population_field, locus_results.hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum);

  return locus_results;

}


// Get a list of hom/het SNPs with a spcified spacing to minimise linkage dis-equilibrium
// and at a specified frequency for the super population. Used as a template for calculating
// the inbreeding coefficient and sample relatedness
std::shared_ptr<const kgl::ContigVariant> kgl::MutationAnalysis::getLocusList( const ContigId_t& contig_id,
                                                                               ContigOffset_t spacing,
                                                                               const std::string& super_pop_freq,
                                                                               double min_frequency,
                                                                               double max_frequency) const {

  // Annotate the variant list with the super population frequency identifier
  std::shared_ptr<ContigVariant> locus_list(std::make_shared<ContigVariant>(super_pop_freq));

  if (unphased_population_->getMap().size() == 1) {

    auto [genome_id, genome_ptr] = *(unphased_population_->getMap().begin());

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

              ExecEnv::log().error("MutationAnalysis::getLocusList, Could not add variant: {}",
                                   variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

            }

          } // for variant array

        } // if spacing

      } // for offset

    } // if contig_opt

  } else {

    ExecEnv::log().error("MutationAnalysis::getLocusList, Unphased Population has {} Genomes, Expected 1",
                         unphased_population_->getMap().size());

  }

  if (locus_list->variantCount() > 0) {

    ExecEnv::log().info("Locus List for super population: {} contains: {} SNPs", locus_list->contigId(), locus_list->variantCount());

  }

  return locus_list;

}

