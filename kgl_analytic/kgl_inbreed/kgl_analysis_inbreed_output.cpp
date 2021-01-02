//
// Created by kellerberrin on 30/11/20.
//


#include "kgl_analysis_inbreed_output.h"
#include "kgl_analysis_inbreed_syngen.h"


#include <fstream>
#include <iostream>


namespace kgl = kellerberrin::genome;


std::string kgl::InbreedingResultColumn::generateIdent(const ContigId_t& contig_id, ContigOffset_t lower, ContigOffset_t upper) {

  std::stringstream ss;

  ss << contig_id << "_" << lower << "_" << upper;

  return ss.str();

}


// Consistency check; all columns have the same row structure.
bool kgl::InbreedParamOutput::verifyResults() const {

  if (column_results_.empty()) {

    return true; // trivally correct.

  }

  // Create a model set of genomes and ensure all the columns have the same structure.
  std::set<GenomeId_t> genome_set;
  auto first_column = column_results_.begin()->results();

  for (auto const& [genome, data] : first_column) {

    genome_set.insert(genome);

  }

  // Check each of the columns for consistency
  for (auto const& column : getColumns()) {

    auto const& column_results = column.results();

    if (genome_set.size() != column_results.size()) {

      ExecEnv::log().warn("InbreedParamOutput::verifyResults, genome column: {} sizes: {} not equal first column size: {}",
                          column.columnIdent(), column_results.size(), genome_set.size());
      return false;

    }

    for (auto const& column : getColumns()) {

      for (auto const& row_result : column.results()) {

        auto const& [row_genome, result] = row_result;
        // If genome_map genome not found then issue a warning and fail.
        auto found_genome = genome_set.find(row_genome);
        if (found_genome == genome_set.end()) {

          ExecEnv::log().warn("InbreedParamOutput::verifyResults, genome ident: {} not found in column: {}", column.columnIdent());
          return false;

        }

      }

    }

  }

  return true;

}


bool kgl::InbreedingOutput::writeNoPedResults(const InbreedParamOutput& output_results, const std::string& file_path) {


  if (output_results.getColumns().empty()) {

    ExecEnv::log().error( "InbreedingOutput::writePedResults; No results to output for param ident",
                          output_results.getParameters().paramIdent());
    return false;

  }

  // Check that the column data is valid.
  if (not output_results.verifyResults()) {

    ExecEnv::log().error( "InbreedingAnalysis::writeColumnResults; Column data corrupt for parameter ident: {}",
                          output_results.getParameters().paramIdent());
    return false;

  }

  // Open the output file.
  std::ofstream outfile;
  std::string output_file_name = Utility::filePath(output_results.getParameters().outputFile(), file_path);
  std::string file_name_ext = output_file_name + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("InbreedingAnalysis::writeColumnResults; could not open output file: {}", file_name_ext);
    return false;

  }


  auto const& parameters = output_results.getParameters();

  outfile << parameters.paramIdent()
          << DELIMITER_ << "Algorithm:" << parameters.inbreedingAlgorthim()
          << DELIMITER_ << "Min_AF:" << parameters.lociiArguments().minAlleleFrequency()
          << DELIMITER_ << "Max_AF:" << parameters.lociiArguments().maxAlleleFrequency()
          << DELIMITER_ << "Spacing:" << parameters.lociiArguments().lociiSpacing()
          << DELIMITER_ << "Count:" << parameters.lociiArguments().lociiCount() << '\n';

  outfile << "Sample";

  for (auto const& column : output_results.getColumns()) {


    outfile << DELIMITER_ << column.columnIdent();

  }

  outfile << '\n';

  std::set<GenomeId_t> genome_set;
  const ResultsMap& first_column_results = output_results.getColumns().front().results();

  for (auto const& [genome, data] : first_column_results) {

    genome_set.insert(genome);

  }

  for (auto const& genome_id : genome_set) {



    outfile << genome_id << DELIMITER_;

    auto columns = output_results.getColumns();
    for (auto const& column : columns) {

      auto find_result = column.results().find(genome_id);
      if (find_result == column.results().end()) {

        ExecEnv::log().error("InbreedingAnalysis::writeColumnResults, Column: {}, Genome sample: {} not found",
                             column.columnIdent(), genome_id);
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


bool kgl::InbreedingOutput::writePedResults( const InbreedParamOutput& output_results,
                                             const GenomePEDData& ped_data,
                                             const std::string& file_path) {

  if (output_results.getColumns().empty()) {

    ExecEnv::log().error("InbreedingOutput::writePedResults; No results to output for parameter ident: {}",
                         output_results.getParameters().paramIdent());
    return false;

  }

  // Check that the column data is valid.
  if (not output_results.verifyResults()) {

    ExecEnv::log().error("InbreedingAnalysis::writeColumnResults; Column data corrupt for parameter ident: {}",
                         output_results.getParameters().paramIdent());
    return false;

  }

  // Open the output file.
  std::ofstream outfile;
  std::string output_file_name = Utility::filePath(output_results.getParameters().outputFile(), file_path);
  std::string file_name_ext = output_file_name + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("InbreedingAnalysis::writeColumnResults; could not open output file: {}", file_name_ext);
    return false;

  }


  auto const& parameters = output_results.getParameters();

  outfile << parameters.paramIdent()
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

  for (auto const& column : output_results.getColumns()) {


    outfile << DELIMITER_ << column.columnIdent();

  }

  outfile << '\n';

  std::set<GenomeId_t> genome_set;
  const ResultsMap& first_column_results = output_results.getColumns().front().results();

  for (auto const& [genome, data] : first_column_results) {

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


    auto columns = output_results.getColumns();
    for (auto const& column : columns) {

      auto find_result = column.results().find(genome_id);
      if (find_result == column.results().end()) {

        ExecEnv::log().error("InbreedingAnalysis::writeColumnResults, Column: {}, Genome sample: {} not found",
                             column.columnIdent(), genome_id);
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



bool kgl::InbreedingOutput::writeSynthetic(const InbreedParamOutput& output_results, const std::string& file_path) {


  if (output_results.getColumns().empty()) {

    ExecEnv::log().error( "InbreedingOutput::writeSynthetic; No results to output for param ident",
                          output_results.getParameters().paramIdent());
    return false;

  }

  // Check that the column data is valid.
  if (not output_results.verifyResults()) {

    ExecEnv::log().error( "InbreedingAnalysis::writeSynthetic; Column data corrupt for parameter ident: {}",
                          output_results.getParameters().paramIdent());
    return false;

  }

  // Open the output file.
  std::ofstream outfile;
  std::string output_file_name = Utility::filePath(output_results.getParameters().outputFile(), file_path);
  std::string file_name_ext = output_file_name + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("InbreedingAnalysis::writeSynthetic; could not open output file: {}", file_name_ext);
    return false;

  }


  auto const& parameters = output_results.getParameters();

  outfile << parameters.paramIdent()
          << DELIMITER_ << "Algorithm:" << parameters.inbreedingAlgorthim()
          << DELIMITER_ << "Min_AF:" << parameters.lociiArguments().minAlleleFrequency()
          << DELIMITER_ << "Max_AF:" << parameters.lociiArguments().maxAlleleFrequency()
          << DELIMITER_ << "Spacing:" << parameters.lociiArguments().lociiSpacing()
          << DELIMITER_ << "Count:" << parameters.lociiArguments().lociiCount() << '\n';

  outfile << "Sample" << DELIMITER_ << "SynInbreed" << DELIMITER_ << "CalcInbreed" << '\n';


  std::set<GenomeId_t> genome_set;
  const ResultsMap& first_column_results = output_results.getColumns().front().results();

  for (auto const& [genome, data] : first_column_results) {

    genome_set.insert(genome);

  }

  for (auto const& genome_id : genome_set) {

    auto [valid, synvalue] = InbreedSynthetic::generateInbreeding(genome_id);

    outfile << genome_id << DELIMITER_ << synvalue << DELIMITER_;

    auto columns = output_results.getColumns();
    for (auto const& column : columns) {

      auto find_result = column.results().find(genome_id);
      if (find_result == column.results().end()) {

        ExecEnv::log().error("InbreedingAnalysis::writeColumnResults, Column: {}, Genome sample: {} not found",
                             column.columnIdent(), genome_id);
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

