//
// Created by kellerberrin on 3/04/23.
//

#include "kgl_pf7_distance_parser.h"
#include "kgl_pf7_sample_parser.h"




namespace kgl = kellerberrin::genome;


double kgl::Pf7DistanceResource::distanceIndex(size_t x_index, size_t y_index) const {

  size_t matrix_index = (y_index * sample_size_) + x_index;

  if (matrix_index >= distance_matrix_size_) {

    ExecEnv::log().warn("Pf7DistanceResource::distanceIndex; x_index: {}, y_index: {}, matrix_size: {}, linear index: {} exceeds matrix size",
                        x_index, y_index, distance_matrix_size_, matrix_index);

    return std::nan(NAN_TEXT_);

  }

  return distance_matrix_ptr_[matrix_index];

}

double kgl::Pf7DistanceResource::distanceIndex(const std::string& sample_x, const std::string& sample_y) const {

  auto x_iter = sample_index_map_ptr_->find(sample_x);
  if (x_iter == sample_index_map_ptr_->end()) {

    ExecEnv::log().warn("Pf7DistanceResource::distanceIndex; x sample id: {} not found", sample_x);
    return std::nan(NAN_TEXT_);

  }

  auto const& [x_id, x_index] = *x_iter;

  auto y_iter = sample_index_map_ptr_->find(sample_y);
  if (y_iter == sample_index_map_ptr_->end()) {

    ExecEnv::log().warn("Pf7DistanceResource::distanceIndex; y sample id: {} not found", sample_y);
    return std::nan(NAN_TEXT_);

  }

  auto const& [y_id, y_index] = *y_iter;

  return distanceIndex(x_index, y_index);

}





bool kgl::ParsePf7Distance::parsePf7Distance(const std::string& matrix_file_name, const std::string& sampleid_file_name) {


  // Get the sample data.
  ParsePf7Sample Pf7_sample_parser;
  if (not Pf7_sample_parser.parsePf7SampleFile(sampleid_file_name)) {

    ExecEnv::log().critical("ParsePf7Distance::parsePf7Distance; failed to create Pf7 Sample Data from file: {}", sampleid_file_name);

  }

  sample_index_map_ptr_->clear();
  size_t index_count{0};
  for (auto const& sample_record : Pf7_sample_parser.getPf7SampleVector()) {

    auto [iter, result] = sample_index_map_ptr_->insert({sample_record.Pf7Sample_id, index_count});
    if (not result) {

      ExecEnv::log().error("ParsePf7Distance::parsePf7Distance; duplicate sample id: {}", sample_record.Pf7Sample_id);

    }
    ++index_count;

  }

  const size_t sample_size = sample_index_map_ptr_->size();
  distance_matrix_size_ = sample_size * sample_size;
  ExecEnv::log().info("ParsePf7Distance::parsePf7Distance; sample id size: {}, distance matrix size: {}",
                      sample_size, distance_matrix_size_);

  ExecEnv::log().info("ParsePf7Distance::parsePf7Distance; reading distance matrix file: {}", matrix_file_name);

  auto distance_line_ptr = SquareTextParser::parseLines(matrix_file_name, distance_matrix_size_);

  ExecEnv::log().info("ParsePf7Distance::parsePf7Distance; completed reading distance matrix file, record count: {}", distance_line_ptr->size());

  if (distance_line_ptr->size() != distance_matrix_size_) {

    ExecEnv::log().critical("ParsePf7Distance::parsePf7Distance; mismatch - distance matrix lines: {}, distance matrix size: {}",
                            distance_line_ptr->size(), distance_matrix_size_);

  }

  ExecEnv::log().info("ParsePf7Distance::parsePf7Distance; begin parsing distance matrix records");

  distance_matrix_ptr_ = std::shared_ptr<double[]>(new double[distance_matrix_size_]);

  for (size_t y_index = 0; y_index < sample_size; ++y_index) {

    for (size_t x_index = 0; x_index < sample_size; ++x_index) {

      size_t vector_index = (y_index * sample_size) + x_index;

      if (vector_index >= distance_matrix_size_) {

        ExecEnv::log().critical("Pf7DistanceResource::distanceIndex; x_index: {}, y_index: {}, matrix_size: {}, linear index: {} exceeds matrix size",
                                x_index, y_index, distance_matrix_size_, vector_index);

      }

      if ((*distance_line_ptr)[vector_index] == NAN_TEXT_) {

        distance_matrix_ptr_[vector_index] = std::nan(NAN_TEXT_);

      } else {

        try {

          double distance_value = std::stod((*distance_line_ptr)[vector_index]);
          if (distance_value < 0) {

            ExecEnv::log().warn("ParsePf7Distance::parsePf7Distance; Distance text: {} not a positive value: {}",
                                (*distance_line_ptr)[vector_index], distance_value);

            distance_matrix_ptr_[vector_index] = std::nan(NAN_TEXT_);

          }

          distance_matrix_ptr_[vector_index] = distance_value;

        } catch(std::exception& e) {

          ExecEnv::log().warn("ParsePf7Distance::parsePf7Distance; Distance text: {} not valid float text, reason: {}, line: {}",
                              (*distance_line_ptr)[vector_index], e.what(), (vector_index + 1));

          distance_matrix_ptr_[vector_index] = std::nan(NAN_TEXT_);

        }

      }

    }

  }

  ExecEnv::log().info("ParsePf7Distance::parsePf7Distance; completed parsing distance matrix records");

  return true;

}


