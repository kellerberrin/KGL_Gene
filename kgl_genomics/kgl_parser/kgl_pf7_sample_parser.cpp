//
// Created by kellerberrin on 31/03/23.
//

#include "kgl_pf7_sample_parser.h"
#include "kel_utility.h"



namespace kgl = kellerberrin::genome;



void kgl::Pf7SampleResource::IndexPf7SampleData() {

  static bool error{false};

  // First load up the canonical records.
  for (auto const& Pf7Sample_record : Pf7Sample_vector_) {

    if (Pf7Sample_record.Pf7Sample_id.empty()) {

      continue;

    }
    auto [iter, result] = Pf7Sample_map_.try_emplace(Pf7Sample_record.Pf7Sample_id, Pf7Sample_record);
    if (not result and not error) {

      ExecEnv::log().warn("Pf7SampleResource::IndexPf7SampleData; duplicate Pf7Sample record ({})",  Pf7Sample_record.Pf7Sample_id);
      error = true;

    }

  }

  ExecEnv::log().info("Pf7SampleResource loaded {}, (Pf7Sample_id, record), lookup pairs", Pf7Sample_map_.size());

}




[[nodiscard]] bool kgl::ParsePf7Sample::parsePf7SampleFile(const std::string& file_name) {


  auto parsed_record_ptr = parseFlatFile(file_name, SquareTextParser::DELIMITER_TSV);

  if (parsed_record_ptr->getRowVector().size() < MINIMUM_ROW_COUNT_) {

    ExecEnv::log().error("ParsePf7Sample::parsePf7SampleFile; Row count: {} for file: {} is below minimum",
                         parsed_record_ptr->getRowVector().size(), file_name);
    return false;

  }

  if (not parsed_record_ptr->checkRowSize(COLUMN_COUNT_)) {

    ExecEnv::log().error("ParsePf7Sample::parsePf7SampleFile; Not all rows have expected column count: {} for file: {}",
                         COLUMN_COUNT_, file_name);
    return false;

  }

  ExecEnv::log().info("Begin Parsing Pf7Sample Gene Resource for file: {}", file_name);

  // Header line is prefixed with '#' and is automatically stripped off.
  size_t record_count{0};
  for (const auto& row_vector :  parsed_record_ptr->getRowVector()) {

    ++record_count;

    // Skip the header.
    if (record_count == 1) {

      continue;

    }

    Pf7SampleRecord Pf7Sample_record;

    Pf7Sample_record.Pf7Sample_id = Utility::trimEndWhiteSpace(row_vector[SAMPLE_ID_OFFSET_]);
    Pf7Sample_record.study_ = Utility::trimEndWhiteSpace(row_vector[STUDY_OFFSET_]);
    Pf7Sample_record.country_ = Utility::trimEndWhiteSpace(row_vector[COUNTRY_OFFSET_]);
    Pf7Sample_record.location1_ = Utility::trimEndWhiteSpace(row_vector[LOCATION1_OFFSET_]);
    Pf7Sample_record.country_latitude_ = Utility::trimEndWhiteSpace(row_vector[COUNTRY_LATITUDE_OFFSET_]);
    Pf7Sample_record.country_longitude_ = Utility::trimEndWhiteSpace(row_vector[COUNTRY_LONGITUDE_OFFSET_]);
    Pf7Sample_record.location1_latitude_ = Utility::trimEndWhiteSpace(row_vector[LOCATION1_LATITUDE_OFFSET_]);
    Pf7Sample_record.location1_longitude_ = Utility::trimEndWhiteSpace(row_vector[LOCATION1_LONGITUDE_OFFSET_]);
    Pf7Sample_record.year_ = Utility::trimEndWhiteSpace(row_vector[YEAR_OFFSET_]);
    Pf7Sample_record.ena_ = Utility::trimEndWhiteSpace(row_vector[ENA_OFFSET_]);
    Pf7Sample_record.all_samples_ = Utility::trimEndWhiteSpace(row_vector[ALL_SAMPLES_OFFSET_]);
    Pf7Sample_record.population_ = Utility::trimEndWhiteSpace(row_vector[POPULATION_OFFSET_]);
    Pf7Sample_record.callable_ = Utility::trimEndWhiteSpace(row_vector[CALLABLE_OFFSET_]);
    Pf7Sample_record.qc_pass_ = Utility::trimEndWhiteSpace(row_vector[QC_PASS_OFFSET_]);
    Pf7Sample_record.qc_fail_reason_ = Utility::trimEndWhiteSpace(row_vector[QC_FAIL_REASON_OFFSET_]);
    Pf7Sample_record.sample_type_ = Utility::trimEndWhiteSpace(row_vector[SAMPLE_TYPE_OFFSET_]);
    Pf7Sample_record.sample_in_pf6_ = Utility::trimEndWhiteSpace(row_vector[SAMPLE_IN_PF6_OFFSET_]);

    Pf7Sample_vector_.push_back(Pf7Sample_record);

  }

  ExecEnv::log().info("ParsePf7Sample::parsePf7SampleFile; Parsed: {} Pf7Sample data records from file: {}",
                      Pf7Sample_vector_.size() , file_name);
  return true;

}
