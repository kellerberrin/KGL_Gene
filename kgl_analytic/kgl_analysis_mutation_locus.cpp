//
// Created by kellerberrin on 13/8/20.
//


#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"


#include <fstream>

namespace kgl = kellerberrin::genome;


// Calculate the HetHom Ratio
bool kgl::InbreedingAnalysis::hetHomRatioLocus() const {

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
    for (auto const&[genome_id, genome_ptr] : diploid_population_->getMap()) {

      auto contig_opt = genome_ptr->getContig(genome_contig_id);

      if (contig_opt) {

        auto result = ped_data_->getMap().find(genome_id);
        if (result == ped_data_->getMap().end()) {

          ExecEnv::log().error("InbreedingAnalysis::hetHomRatioLocu, Genome sample: {} does not have a PED record", genome_id);
          continue;

        }
        auto const& [sample_id, ped_record] = *result;
        auto locus_result = locus_map.find(ped_record.superPopulation());
        if (locus_result == locus_map.end()) {

          ExecEnv::log().error("InbreedingAnalysis::hetHomRatioLocu, Locus set not found for super population: {}", ped_record.superPopulation());
          continue;

        }
        auto const& [super_pop_id, locus_list] = *locus_result;

        std::future<LocusResults> future = thread_pool.enqueueTask(&InbreedingAnalysis::multiLocus1,
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
bool kgl::InbreedingAnalysis::syntheticInbreeding() const {

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


          std::future<LocusResults> future = thread_pool.enqueueTask(&InbreedingAnalysis::processRitlandLocus,
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



bool kgl::InbreedingAnalysis::syntheticResults( const ContigId_t& contig_id,
                                              const std::map<GenomeId_t, LocusResults>& genome_results_map) const {

  // Append the results.
  std::ofstream outfile;
  outfile.open(output_file_name_, std::ofstream::out |  std::ofstream::app);

  outfile << contig_id << DELIMITER_ << "SuperPop" << DELIMITER_ << "Inbreeding"
          << DELIMITER_ << "HetCount" << DELIMITER_ << "HomCount"
          << DELIMITER_ << "Het/Hom"<< DELIMITER_ << "TotalLoci" << DELIMITER_ << "Ritland\n";

  for (auto const& [genome_id, locus_results] : genome_results_map) {


    auto const [valid_value, inbreeding] = generateInbreeding(genome_id);
    double het_hom_ratio = (locus_results.homo_count > 0 ? static_cast<double>(locus_results.hetero_count)/static_cast<double>(locus_results.homo_count) : 0.0);

    outfile << genome_id << DELIMITER_;
    outfile << genome_id.substr(0, genome_id.find_first_of("_")) << DELIMITER_;

    if (valid_value) {

      outfile << inbreeding << DELIMITER_;

    } else {

      outfile << "Invalid" << DELIMITER_;

    }
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


bool kgl::InbreedingAnalysis::writeResults( const ContigId_t& contig_id,
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

      ExecEnv::log().error("InbreedingAnalysis::hetHomRatio, Genome sample: {} does not have a PED record", genome_id);
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



// Get a list of hom/het SNPs with a spcified spacing to minimise linkage dis-equilibrium
// and at a specified frequency for the super population. Used as a template for calculating
// the inbreeding coefficient and sample relatedness
std::shared_ptr<const kgl::ContigVariant> kgl::InbreedingAnalysis::getLocusList( const ContigId_t& contig_id,
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

              ExecEnv::log().error("InbreedingAnalysis::getLocusList, Could not add variant: {}",
                                   variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

            }

          } // for variant array

        } // if spacing

      } // for offset

    } // if contig_opt

  } else {

    ExecEnv::log().error("InbreedingAnalysis::getLocusList, Unphased Population has {} Genomes, Expected 1",
                         unphased_population_->getMap().size());

  }

  if (locus_list->variantCount() > 0) {

    ExecEnv::log().info("Locus List for super population: {} contains: {} SNPs", locus_list->contigId(), locus_list->variantCount());

  }

  return locus_list;

}

