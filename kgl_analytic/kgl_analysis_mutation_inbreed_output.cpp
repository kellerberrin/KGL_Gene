//
// Created by kellerberrin on 30/11/20.
//


#include "kgl_analysis_mutation_inbreed_output.h"
#include "kgl_analysis_mutation_syngen.h"


#include <fstream>
#include <iostream>
#include <iomanip>


namespace kgl = kellerberrin::genome;




bool kgl::InbreedingOutputResults::verifyResults() const {

  if (results_vector_.empty()) {

    return true; // trivally correct.

  }

  // Create a model set of genomes and ensure all the columns have the same structure.
  std::set<GenomeId_t> genome_set;
  auto [ident, first_column] = *results_vector_.begin();

  for (auto const& [genome, data] : first_column) {

    genome_set.insert(genome);

  }

  // Check each of the columns for consistency
  for (auto const& [parameters, genome_map] : results_vector_) {

    if (genome_set.size() != genome_map.size()) {

      ExecEnv::log().warn("InbreedingAnalysis::verifyColumnMap, genome column sizes: {} not equal first column size: {}",
                          genome_map.size(), genome_set.size());
      return false;

    }

    for (auto const& [col_genome, col_data] : genome_map) {

      // If genome_map genome not found then issue a warning and fail.
      auto result = genome_set.find(col_genome);
      if (result == genome_set.end()) {

        ExecEnv::log().warn("InbreedingAnalysis::verifyColumnMap, genome ident: {} not found in first column", col_genome);
        return false;

      }

    }

  }

  return true;

}



bool kgl::InbreedingOutput::writeSynResults( const InbreedingOutputResults& column_results,
                                             const std::string& file_path) {


  if (column_results.resultsVector().empty()) {

    ExecEnv::log().error("InbreedingOutput::writeSynResults; results vector empty for identifier: {}", column_results.identifier());
    return false;

  }

  // Check that the column data is valid.
  if (not column_results.verifyResults()) {

    ExecEnv::log().error("InbreedingOutput::writeSynResults; Column data corrupt in result identifier: {}", column_results.identifier());
    return false;

  }

  // Open the output file.
  std::ofstream outfile;
  std::string file_name_ext = file_path + column_results.identifier() + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out | std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("IInbreedingOutput::writeSynResults; could not open output file: {}", file_name_ext);
    return false;

  }

  // Sort the column results using map
  std::map<std::string, ResultsMap> column_map;
  for (auto const&[parameters, results] : column_results.resultsVector()) {

    std::stringstream ss;

    ss << std::setw(9) << std::setfill('0')
       << parameters.lociiArguments().lowerOffset()
       << "_" << std::setw(9) << std::setfill('0')
       << parameters.lociiArguments().upperOffset();

    column_map[ss.str()] = results;

  }

  auto[parameters, results] = column_results.resultsVector().front();

  outfile << column_results.identifier()
          << DELIMITER_ << "Algorithm:" << parameters.inbreedingAlgorthim()
          << DELIMITER_ << "Min_AF:" << parameters.lociiArguments().minAlleleFrequency()
          << DELIMITER_ << "Max_AF:" << parameters.lociiArguments().maxAlleleFrequency()
          << DELIMITER_ << "Spacing:" << parameters.lociiArguments().lociiSpacing()
          << DELIMITER_ << "Count:" << parameters.lociiArguments().lociiCount() << '\n';

  outfile << "Sample" << DELIMITER_
          << "SynInbreed";

  for (auto const&[column_id, genome_results] : column_map) {


    outfile << DELIMITER_ << column_id;

  }

  outfile << '\n';

  std::set<GenomeId_t> genome_set;
  auto[ident, first_column] = *column_map.begin();

  for (auto const&[genome, data] : first_column) {

    genome_set.insert(genome);

  }

  for (auto const &genome_id : genome_set) {

    outfile << genome_id << DELIMITER_;

    auto [valid_flag, syn_inbreed] = InbreedSynthetic::generateInbreeding(genome_id);

    outfile << syn_inbreed;

    for (auto const&[column_id, genome_results] : column_map) {

      auto find_result = genome_results.find(genome_id);
      if (find_result == genome_results.end()) {

        ExecEnv::log().error("InbreedingAnalysis::writeColumnResults, Column: {}, Genome sample: {} not found", column_id, genome_id);
        return false;

      }

      auto[genome_id, inbreed_data] = *find_result;
      outfile << DELIMITER_ << inbreed_data.inbred_allele_sum ;

    }

    outfile << '\n';

  }

  outfile.flush();

  return true;

}


bool kgl::InbreedingOutput::writeColumnResults( const InbreedingOutputResults& column_results,
                                                const GenomePEDData& ped_data,
                                                const std::string& file_path) {

  if (column_results.resultsVector().empty()) {

    ExecEnv::log().error("InbreedingAnalysis::writeColumnResults; results vector empty for identifier: {}", column_results.identifier());
    return false;

  }

  // Check that the column data is valid.
  if (not column_results.verifyResults()) {

    ExecEnv::log().error("InbreedingAnalysis::writeColumnResults; Column data corrupt in result identifier: {}", column_results.identifier());
    return false;

  }

  // Open the output file.
  std::ofstream outfile;
  std::string file_name_ext = file_path + column_results.identifier() + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("InbreedingAnalysis::writeColumnResults; could not open output file: {}", file_name_ext);
    return false;

  }

  // Sort the column results using map
  std::map<std::string, ResultsMap> column_map;
  for (auto const& [parameters, results] : column_results.resultsVector()) {

    std::stringstream ss;

    ss << std::setw(9) << std::setfill('0')
       << parameters.lociiArguments().lowerOffset()
       << "_" << std::setw(9) << std::setfill('0')
       << parameters.lociiArguments().upperOffset();

    column_map[ss.str()] = results;

  }

  auto [parameters, results] = column_results.resultsVector().front();

  outfile << column_results.identifier()
          << DELIMITER_ << "Algorithm:" << parameters.inbreedingAlgorthim()
          << DELIMITER_ << "Min_AF:" << parameters.lociiArguments().minAlleleFrequency()
          << DELIMITER_ << "Max_AF:" << parameters.lociiArguments().maxAlleleFrequency()
          << DELIMITER_ << "Spacing:" << parameters.lociiArguments().lociiSpacing()
          << DELIMITER_ << "Count:" << parameters.lociiArguments().lociiCount() << '\n';

  outfile << "Sample"
          << DELIMITER_ << "Population"
          << DELIMITER_ << "Description"
          << DELIMITER_ << "SuperPopulation"
          << DELIMITER_ << "Description"
          << DELIMITER_ << "Relationship"
          << DELIMITER_ << "Sex"
          << DELIMITER_ << "Mother"
          << DELIMITER_ << "Father";

  for (auto const& [column_id, genome_results] : column_map) {


    outfile << DELIMITER_ << column_id;

  }

  outfile << '\n';

  std::set<GenomeId_t> genome_set;
  auto [ident, first_column] = *column_map.begin();

  for (auto const& [genome, data] : first_column) {

    genome_set.insert(genome);

  }

  for (auto const& genome_id : genome_set) {

    auto result = ped_data.getMap().find(genome_id);

    if (result == ped_data.getMap().end()) {

      ExecEnv::log().error("InbreedingAnalysis::writeColumnResults, Genome sample: {} does not have a PED record", genome_id);
      continue;

    }

    auto const& [sample_id, ped_record] = *result;

    outfile << genome_id << DELIMITER_;
    outfile << ped_record.population() << DELIMITER_;
    outfile << ped_record.populationDescription() << DELIMITER_;
    outfile << ped_record.superPopulation() << DELIMITER_;
    outfile << ped_record.superDescription() << DELIMITER_;
    outfile << ped_record.relationship() << DELIMITER_;
    outfile << ped_record.sex() << DELIMITER_;
    outfile << ped_record.maternalId() << DELIMITER_;
    outfile << ped_record.paternalId() << DELIMITER_;


    for (auto const& [column_id, genome_results] : column_map) {

      auto find_result = genome_results.find(genome_id);
      if (find_result == genome_results.end()) {

        ExecEnv::log().error("InbreedingAnalysis::writeColumnResults, Column: {}, Genome sample: {} not found", column_id, genome_id);
        return false;

      }

      auto [genome_id, inbreed_data] = *find_result;
      outfile << inbreed_data.inbred_allele_sum << DELIMITER_;

    }

    outfile << '\n';

  }

  outfile.flush();

  return true;

}


bool kgl::InbreedingOutput::writeResults( const ResultsMap& genome_results_map,
                                            const GenomePEDData& ped_data,
                                            const std::string& output_file_name,
                                            InbreedingParameters& parameters) {


  // Open the output file.
  std::ofstream outfile;
  std::string file_name_ext = output_file_name + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreedingSample; could not open output file: {}", file_name_ext);

  }

  outfile << "Min_AF:" << parameters.lociiArguments().minAlleleFrequency()
          << DELIMITER_ << "Max_AF:" << parameters.lociiArguments().maxAlleleFrequency()
          << DELIMITER_ << "Spacing:" << parameters.lociiArguments().lociiSpacing()
          << DELIMITER_ << "LowerOffset:" << parameters.lociiArguments().lowerOffset()
          << DELIMITER_ << "UpperOffset:" << parameters.lociiArguments().upperOffset() << '\n';

  outfile << DELIMITER_ << "Population"
          << DELIMITER_ << "Description"
          << DELIMITER_ << "SuperPopulation"
          << DELIMITER_ << "Description"
          << DELIMITER_ << "MajorHet"
          << DELIMITER_ << "MajorHetFreq"
          << DELIMITER_ << "MinorHet"
          << DELIMITER_ << "MinorHetFreq"
          << DELIMITER_ << "MajorHom"
          << DELIMITER_ << "MajorHomFreq"
          << DELIMITER_ << "MinorHom"
          << DELIMITER_ << "MinorHomFreq"
          << DELIMITER_ << "Het/Hom"
          << DELIMITER_ << "TotalLoci"
          << DELIMITER_ << "CalcInbreed" << '\n';

  for (auto const& [genome_id, locus_results] : genome_results_map) {

    auto result = ped_data.getMap().find(genome_id);

    if (result == ped_data.getMap().end()) {

      ExecEnv::log().error("InbreedingAnalysis::writeResults, Genome sample: {} does not have a PED record", genome_id);
      continue;

    }

    auto const& [sample_id, ped_record] = *result;

    double het_hom_ratio;
    if (locus_results.minor_homo_count > 0) {

      size_t total_het = locus_results.major_hetero_count + locus_results.minor_hetero_count;
      het_hom_ratio = static_cast<double>(total_het) / static_cast<double>(locus_results.minor_homo_count);

    } else {

      het_hom_ratio = 0.0;

    }

    outfile << genome_id << DELIMITER_;
    outfile << ped_record.population() << DELIMITER_;
    outfile << ped_record.populationDescription() << DELIMITER_;
    outfile << ped_record.superPopulation() << DELIMITER_;
    outfile << ped_record.superDescription() << DELIMITER_;
    outfile << locus_results.major_hetero_count << DELIMITER_;
    outfile << locus_results.major_hetero_freq << DELIMITER_;
    outfile << locus_results.minor_hetero_count << DELIMITER_;
    outfile << locus_results.minor_hetero_freq << DELIMITER_;
    outfile << locus_results.major_homo_count << DELIMITER_;
    outfile << locus_results.major_homo_freq << DELIMITER_;
    outfile << locus_results.minor_homo_count << DELIMITER_;
    outfile << locus_results.minor_homo_freq << DELIMITER_;
    outfile << het_hom_ratio << DELIMITER_;
    outfile << locus_results.total_allele_count << DELIMITER_;
    outfile << locus_results.inbred_allele_sum << DELIMITER_;
    outfile << '\n';

  }

  outfile.flush();

  return true;

}

