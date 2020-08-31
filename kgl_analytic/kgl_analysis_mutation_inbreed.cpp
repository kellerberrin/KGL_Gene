//
// Created by kellerberrin on 13/8/20.
//


#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_analysis_mutation_inbreed_aux.h"


#include <fstream>

namespace kgl = kellerberrin::genome;


// Calculate the Population Inbreeding Coefficient
bool kgl::InbreedingAnalysis::populationInbreeding(std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                   const DiploidPopulation& diploid_population,
                                                   const GenomePEDData& ped_data,
                                                   const std::string& output_file_name) {

  // Open the output file.
  std::ofstream outfile5;
  std::string file_name_ext = output_file_name + "pop_All" + FILE_EXT_;
  outfile5.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile5 << diploid_population.populationId() << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.0
           << DELIMITER_ << "Max_AF:" << 1.00 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                                           0.0,   // min af
                                                                           1.0,   // max af
                                                                           1000);  // spacing

  processResults(contig_locus_map, diploid_population, ped_data, outfile5);

  // Open the output file.
  std::ofstream outfile1;
  file_name_ext = output_file_name + "pop_01" + FILE_EXT_;
  outfile1.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile1 << diploid_population.populationId() << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.0
           << DELIMITER_ << "Max_AF:" << 0.01 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                            0.0,   // min af
                                                            0.01,   // max af
                                                            1000);  // spacing

  processResults(contig_locus_map, diploid_population, ped_data, outfile1);

  // Open the output file.
  std::ofstream outfile2;
  file_name_ext = output_file_name + "pop_05" + FILE_EXT_;
  outfile2.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile2 << diploid_population.populationId() << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.01
           << DELIMITER_ << "Max_AF:" << 0.05 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                            0.01,   // min af
                                                            0.05,   // max af
                                                            1000);  // spacing

  processResults(contig_locus_map, diploid_population, ped_data, outfile2);

  // Open the output file.
  std::ofstream outfile3;
  file_name_ext = output_file_name + "pop_20" + FILE_EXT_;
  outfile3.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile3 << diploid_population.populationId() << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.05
           << DELIMITER_ << "Max_AF:" << 0.20 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                            0.05,   // min af
                                                            0.2,   // max af
                                                            1000);  // spacing

  processResults(contig_locus_map, diploid_population, ped_data, outfile3);

  // Open the output file.
  std::ofstream outfile4;
  file_name_ext = output_file_name + "pop_100" + FILE_EXT_;
  outfile4.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile4 << diploid_population.populationId() << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.20
           << DELIMITER_ << "Max_AF:" << 1.00 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                            0.2,   // min af
                                                            1.0,   // max af
                                                            1000);  // spacing

  processResults(contig_locus_map, diploid_population, ped_data, outfile4);



  return true;

}


bool kgl::InbreedingAnalysis::processResults( const ContigLocusMap& contig_locus_map,
                                              const DiploidPopulation& diploid_population,
                                              const GenomePEDData& ped_data,
                                              std::ostream& outfile) {

  for (auto const& [genome_contig_id, locus_map] : contig_locus_map) {

    // Use a thread pool to calculate inbreeding and relatedness.
    ThreadPool thread_pool;
    std::vector<std::future<LocusResults>> future_vector;
    for (auto const&[genome_id, genome_ptr] : diploid_population.getMap()) {

      auto contig_opt = genome_ptr->getContig(genome_contig_id);

      if (contig_opt) {

        auto result = ped_data.getMap().find(genome_id);
        if (result == ped_data.getMap().end()) {

          ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Genome sample: {} does not have a PED record", genome_id);
          continue;

        }
        auto const& [sample_id, ped_record] = *result;
        auto locus_result = locus_map.find(ped_record.superPopulation());
        if (locus_result == locus_map.end()) {

          ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Locus set not found for super population: {}", ped_record.superPopulation());
          continue;

        }
        auto const& [super_pop_id, locus_list] = *locus_result;

        std::future<LocusResults> future = thread_pool.enqueueTask(&InbreedingCalculation::processRitlandLocus,
                                                                   genome_id,
                                                                   contig_opt.value(),
                                                                   InbreedSampling::lookupSuperPopulationField(super_pop_id),
                                                                   locus_list);
        future_vector.push_back(std::move(future));

      }

    }

    // Retrieve the thread results into a map.
    ResultsMap genome_results_map;
    for (auto &future : future_vector) {

      auto locus_results = future.get();
      ExecEnv::log().info("Processed Diploid genome: {} for inbreeding and relatedness", locus_results.genome);
      genome_results_map[locus_results.genome] = locus_results;


    }

    if (not genome_results_map.empty()) {

      writeResults(genome_contig_id, genome_results_map, ped_data, outfile);

    }

  }

  return true;

}


bool kgl::InbreedingAnalysis::writeResults( const ContigId_t& contig_id,
                                            const ResultsMap& genome_results_map,
                                            const GenomePEDData& ped_data,
                                            std::ostream& outfile) {


  outfile << contig_id << DELIMITER_ << "Population" << DELIMITER_
          << "Description" << DELIMITER_ << "SuperPopulation" << DELIMITER_ << "Description" << DELIMITER_
          << "HetCount" << DELIMITER_ << "HomCount" << DELIMITER_ << "Het/Hom"<< DELIMITER_
          << "TotalLoci" << DELIMITER_ << "Ritland\n";

  for (auto const& [genome_id, locus_results] : genome_results_map) {

    auto result = ped_data.getMap().find(genome_id);

    if (result == ped_data.getMap().end()) {

      ExecEnv::log().error("InbreedingAnalysis::writeResults, Genome sample: {} does not have a PED record", genome_id);
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


// Calculate the Synthetic Population Inbreeding Coefficient
bool kgl::InbreedingAnalysis::syntheticInbreeding(  std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                    const std::string& output_file_name) {

  // Open the output file.
  std::ofstream outfile5;
  std::string file_name_ext = output_file_name + "syn_All" + FILE_EXT_;
  outfile5.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile5 << "Synthetic" << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.00
           << DELIMITER_ << "Max_AF:" << 1.00 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                                           0.0,   // min af
                                                                           1.0,   // max af
                                                                           1000);  // spacing

  processSynResults(contig_locus_map, outfile5);

  // Open the output file.
  std::ofstream outfile1;
  file_name_ext = output_file_name + "syn_01" + FILE_EXT_;
  outfile1.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile1 << "Synthetic" << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.0
          << DELIMITER_ << "Max_AF:" << 0.01 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                            0.0,   // min af
                                                            0.01,   // max af
                                                            1000);  // spacing

  processSynResults(contig_locus_map, outfile1);

  // Open the output file.
  std::ofstream outfile2;
  file_name_ext = output_file_name + "syn_05" + FILE_EXT_;
  outfile2.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile2 << "Synthetic" << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.01
          << DELIMITER_ << "Max_AF:" << 0.05 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                            0.01,   // min af
                                                            0.05,   // max af
                                                            1000);  // spacing

  processSynResults(contig_locus_map, outfile2);

  // Open the output file.
  std::ofstream outfile3;
  file_name_ext = output_file_name + "syn_20" + FILE_EXT_;
  outfile3.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile3 << "Synthetic" << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.05
           << DELIMITER_ << "Max_AF:" << 0.20 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                            0.05,   // min af
                                                            0.2,   // max af
                                                            1000);  // spacing

  processSynResults(contig_locus_map, outfile3);

  // Open the output file.
  std::ofstream outfile4;
  file_name_ext = output_file_name + "syn_100" + FILE_EXT_;
  outfile4.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile4 << "Synthetic" << DELIMITER_ << "Min_AF:" << DELIMITER_ << 0.20
           << DELIMITER_ << "Max_AF:" << 1.00 << DELIMITER_ << "Spacing:" << 1000 << '\n';

  contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                            0.2,   // min af
                                                            1.0,   // max af
                                                            1000);  // spacing

  processSynResults(contig_locus_map, outfile4);


  return true;

}


bool kgl::InbreedingAnalysis::processSynResults( const ContigLocusMap& contig_locus_map,
                                                 std::ostream& outfile) {

  // For each contig.
  for (auto const& [genome_contig_id, locus_map] : contig_locus_map) {

    // Use a thread pool to calculate inbreeding and relatedness.
    for (auto const&[super_pop_id, locus_list] : locus_map) {

      std::shared_ptr<const DiploidPopulation> population = InbreedSampling::generateSyntheticPopulation( MIN_INBREEDING_COEFICIENT,
                                                                                                          MAX_INBREEDING_COEFICIENT,
                                                                                                          STEP_INBREEDING_COEFICIENT,
                                                                                                          super_pop_id,
                                                                                                          *locus_list);

      ThreadPool thread_pool;
      std::vector<std::future<LocusResults>> future_vector;
      for (auto const&[genome_id, genome_ptr] : population->getMap()) {

        auto contig_opt = genome_ptr->getContig(genome_contig_id);

        if (contig_opt) {


          std::future<LocusResults> future = thread_pool.enqueueTask(&InbreedingCalculation::processRitlandLocus,
                                                                     genome_id,
                                                                     contig_opt.value(),
                                                                     InbreedSampling::lookupSuperPopulationField(super_pop_id),
                                                                     locus_list);
          future_vector.push_back(std::move(future));

        }

      }

      // Retrieve the thread results into a map.
      ResultsMap genome_results_map;
      for (auto &future : future_vector) {

        auto locus_results = future.get();
        ExecEnv::log().info("Processed Diploid genome: {} for inbreeding and relatedness", locus_results.genome);
        genome_results_map[locus_results.genome] = locus_results;


      }

      if (not genome_results_map.empty()) {

        writeSynResults(genome_contig_id, genome_results_map, outfile);

      }

    } // for locus

  } // for contig.

  return true;

}


bool kgl::InbreedingAnalysis::writeSynResults(const ContigId_t& contig_id,
                                              const std::map<GenomeId_t, LocusResults>& genome_results_map,
                                              std::ostream& outfile) {


  outfile << contig_id << DELIMITER_ << "SuperPop" << DELIMITER_ << "Inbreeding"
          << DELIMITER_ << "HetCount" << DELIMITER_ << "HomCount"
          << DELIMITER_ << "Het/Hom"<< DELIMITER_ << "TotalLoci" << DELIMITER_ << "Ritland\n";

  for (auto const& [genome_id, locus_results] : genome_results_map) {


    auto const [valid_value, inbreeding] = InbreedSampling::generateInbreeding(genome_id);
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

