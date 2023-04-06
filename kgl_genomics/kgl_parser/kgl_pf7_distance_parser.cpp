//
// Created by kellerberrin on 3/04/23.
//

#include "kgl_pf7_distance_parser.h"
#include "kgl_pf7_sample_parser.h"
#include "kel_mt_buffer.h"

#include <limits>

namespace kgl = kellerberrin::genome;


double kgl::Pf7DistanceResource::getDistance(size_t x_index, size_t y_index) const {


  if (x_index >= sample_size_ or y_index >= sample_size_) {

    ExecEnv::log().warn("Pf7DistanceResource::getDistance; x_index: {}, y_index: {} out of bounds, indices must be in the range [0, {})",
                        x_index, y_index, sample_size_);

    return std::nan(NAN_TEXT_);

  }

  size_t matrix_index = (y_index * sample_size_) + x_index;

  return distance_matrix_ptr_[matrix_index];

}

double kgl::Pf7DistanceResource::getDistance(const std::string& sample_x, const std::string& sample_y) const {

  auto x_iter = sample_index_map_ptr_->find(sample_x);
  if (x_iter == sample_index_map_ptr_->end()) {

    ExecEnv::log().warn("Pf7DistanceResource::getDistance; x sample id: {} not found", sample_x);
    return std::nan(NAN_TEXT_);

  }

  auto const& [x_id, x_index] = *x_iter;

  auto y_iter = sample_index_map_ptr_->find(sample_y);
  if (y_iter == sample_index_map_ptr_->end()) {

    ExecEnv::log().warn("Pf7DistanceResource::getDistance; y sample id: {} not found", sample_y);
    return std::nan(NAN_TEXT_);

  }

  auto const& [y_id, y_index] = *y_iter;

  return getDistance(x_index, y_index);

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

  // Allocate the distance matrix
  auto matrix_ptr = new (std::nothrow) double [distance_matrix_size_];
  if (matrix_ptr == nullptr) {

    ExecEnv::log().critical("ParsePf7Distance::parsePf7Distance; failed to memory allocate matrix size: {}", distance_matrix_size_);

  }
  distance_matrix_ptr_ = std::shared_ptr<double[]>(matrix_ptr);

  // Open the distance matrix file.
  auto file_io_opt = BaseStreamIO::getStreamIO(matrix_file_name);
  if (not file_io_opt) {

    ExecEnv::log().critical("ParsePf7Distance::parsePf7Distance; I/O error; could not open file: {}", matrix_file_name);

  }
  auto file_io_ptr = std::move(file_io_opt.value());

  ExecEnv::log().info("ParsePf7Distance::parsePf7Distance; begin parsing distance matrix records");

  // Read each record and place in the appropriate matrix element.
  size_t vector_index{0};
  while (true) {

    auto line_record = file_io_ptr->readLine();
    if (line_record.EOFRecord()) break;

    auto [line_number, line_string] = line_record.getLineData();

    if (vector_index >= distance_matrix_size_) {

      ExecEnv::log().critical("ParsePf7Distance::parsePf7Distance; matrix_size: {}, linear index: {} exceeds matrix size",
                              distance_matrix_size_, vector_index);

    }

    if (line_string == NAN_TEXT_) {

      distance_matrix_ptr_[vector_index] = std::nan(NAN_TEXT_);

    } else {

      try {

        double distance_value = std::stod(line_string);
        if (distance_value < 0) {

          ExecEnv::log().warn("ParsePf7Distance::parsePf7Distance; Distance text: {} not a positive value: {}",
                              line_string, distance_value);

          distance_matrix_ptr_[vector_index] = std::nan(NAN_TEXT_);

        }

        distance_matrix_ptr_[vector_index] = distance_value;

      } catch(std::exception& e) {

        ExecEnv::log().warn("ParsePf7Distance::parsePf7Distance; Distance text: {} not valid float text, reason: {}, line: {}",
                            line_string, e.what(), (vector_index + 1));

        distance_matrix_ptr_[vector_index] = std::nan(NAN_TEXT_);

      }

    }

    ++vector_index;

  }

  ExecEnv::log().info("ParsePf7Distance::parsePf7Distance; completed parsing distance matrix records");

  return true;

}


