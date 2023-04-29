//
// Created by kellerberrin on 3/04/23.
//

#include "kgl_pf7_fws_parser.h"


namespace kgl = kellerberrin::genome;



void kgl::Pf7FwsResource::indexPf7FwsData() {

  // First load up the canonical records.
  for (auto const& Pf7Fws_record : Pf7Fws_vector_) {

    if (Pf7Fws_record.Pf7Sample_id.empty()) {

      continue;

    }
    auto [iter, result] = Pf7Fws_map_.try_emplace(Pf7Fws_record.Pf7Sample_id, Pf7Fws_record);
    if (not result) {

      ExecEnv::log().warn("Pf7SampleResource::IndexPf7SampleData; duplicate Pf7Sample record ({})",  Pf7Fws_record.Pf7Sample_id);

    }

  }

  ExecEnv::log().info("Pf7FwsResource loaded {}, (Pf7Sample_id, record), lookup pairs", Pf7Fws_map_.size());

}

double kgl::Pf7FwsResource::getFWS(const GenomeId_t& genome_id) const {

  if (getMap().contains(genome_id)) {

    auto sample_iter = getMap().find(genome_id);
    auto const &[id, sample_data] = *sample_iter;
    return sample_data.FWS_value;

  }

  ExecEnv::log().warn("Pf7FwsResource::filterFWS; Unable to find FWS statistic for genome: {}", genome_id);
  return std::nan("n/a");

}


// Population must be PF7
// Important - this code above only filters a shallow copy of the population.
std::shared_ptr<kgl::PopulationDB> kgl::Pf7FwsResource::filterFWS( FwsFilterType filter_type,
                                                                   double fws_threshold,
                                                                   const std::shared_ptr<const PopulationDB>& Pf7_unfiltered_ptr) const {

  auto filtered_ptr = std::make_shared<kgl::PopulationDB>(Pf7_unfiltered_ptr->populationId() + "_FWS_Filtered",
                                                          Pf7_unfiltered_ptr->dataSource());

  for (auto const& [genome_id, genome_ptr] : Pf7_unfiltered_ptr->getMap()) {

    if (getMap().contains(genome_id)) {

      auto sample_iter = getMap().find(genome_id);
      auto const& [id, sample_data] = *sample_iter;
      bool accept = filter_type == FwsFilterType::GREATER_EQUAL ? sample_data.FWS_value >= fws_threshold : sample_data.FWS_value <= fws_threshold;
      if (accept) {

        if (not filtered_ptr->addGenome(genome_ptr)) {

          ExecEnv::log().warn("Pf7FwsResource::filterFWS; Unable to add filtered genome: {} to filtered  population", genome_id);

        }

      } // Pass

    } else {

      ExecEnv::log().warn("PPf7FwsResource::filterFWS; Genome: {} not found in sample data", genome_id);

    }

  } // For genomes.

  return filtered_ptr;

}


[[nodiscard]] bool kgl::ParsePf7Fws::parsePf7FwsFile(const std::string& file_name) {


  auto parsed_record_ptr = parseFlatFile(file_name, SquareTextParser::DELIMITER_TSV);

  if (parsed_record_ptr->getRowVector().size() < MINIMUM_ROW_COUNT_) {

    ExecEnv::log().error("ParsePf7Fws::parsePf7FwsFile; Row count: {} for file: {} is below minimum",
                         parsed_record_ptr->getRowVector().size(), file_name);
    return false;

  }

  if (not parsed_record_ptr->checkRowSize(COLUMN_COUNT_)) {

    ExecEnv::log().error("ParsePf7Fws::parsePf7FwsFile; Not all rows have expected column count: {} for file: {}",
                         COLUMN_COUNT_, file_name);
    return false;

  }

  ExecEnv::log().info("Begin Parsing Pf7 within-host infection fixation index (FWS) resource for file: {}", file_name);

  size_t record_count{0};
  for (const auto& row_vector :  parsed_record_ptr->getRowVector()) {

    ++record_count;

    // Skip the header.
    if (record_count == 1) {

      continue;

    }

    Pf7FwsRecord Pf7Fws_record;

    Pf7Fws_record.Pf7Sample_id = Utility::trimEndWhiteSpace(row_vector[SAMPLE_ID_OFFSET_]);
    auto FWS_float_text = Utility::trimEndWhiteSpace(row_vector[FWS_OFFSET_]);

    try {

      double FWS_Value = std::stod(FWS_float_text);
      Pf7Fws_record.FWS_value = FWS_Value;

    } catch(std::exception& e) {

      ExecEnv::log().info("ParsePf7Fws::parsePf7FwsFile; FWS text: {} not valid float text, reason: {}, line: {}, file: {}",
                          FWS_float_text, e.what(), record_count, file_name);
      continue;

    }

    Pf7Fws_vector_.push_back(Pf7Fws_record);

  }

  ExecEnv::log().info("ParsePf7Fws::parsePf7FwsFile; Parsed: {} Pf7 FWS data records from file: {}",
                      Pf7Fws_vector_.size() , file_name);
  return true;

}


