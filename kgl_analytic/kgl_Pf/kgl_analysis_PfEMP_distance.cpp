//
// Created by kellerberrin on 6/04/23.
//


#include "kgl_analysis_PfEMP.h"


namespace kgl = kellerberrin::genome;



void kgl::PfEMPAnalysis::checkDistanceMatrix( const std::shared_ptr<const PopulationDB>& all_population_ptr,
                                              const std::shared_ptr<const PopulationDB>& filtered_population_ptr) const {

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Diagonal distance
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  const size_t distance_sample_size = Pf7_distance_ptr_->distanceSamples()->size();
  double sum_diagonal{0.0};
  size_t diagonal_count{0};
  size_t nan_count{0};
  for (size_t x_index = 0; x_index < distance_sample_size; ++x_index) {

    double distance = Pf7_distance_ptr_->getDistance(x_index, x_index);
    if (std::isnan(distance)) {

      ++nan_count;

    } else {

      sum_diagonal += distance;
      ++diagonal_count;

    }

  }

  double av_distance = sum_diagonal / static_cast<double>(diagonal_count);
  ExecEnv::log().info("Using Index; Diagonal average distance: {}, Diagonal non-nan: {}, Diagonal nan: {}", av_distance, diagonal_count, nan_count);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Examine sample id lookup.
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  sum_diagonal = 0.0;
  diagonal_count = 0;
  nan_count = 0;
  for (auto const& [genome_id, genome_ptr] : all_population_ptr->getMap()) {

    double distance = Pf7_distance_ptr_->getDistance(genome_id, genome_id);
    if (std::isnan(distance)) {

      ++nan_count;

    } else {

      sum_diagonal += distance;
      ++diagonal_count;

    }

  }
  av_distance = sum_diagonal / static_cast<double>(diagonal_count);
  ExecEnv::log().info("Using Sample Id; Diagonal average distance: {}, Diagonal non-nan: {}, Diagonal nan: {}", av_distance, diagonal_count, nan_count);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Check the smallest each row is self.
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::pair<double, std::string> smallest{ 1000.0, ""};
  size_t self_smallest{0};
  size_t other_smallest{0};

  for (auto const& [y_genome_id, y_genome_ptr] : filtered_population_ptr->getMap()) {

    for (auto const& [x_genome_id, x_genome_ptr] : filtered_population_ptr->getMap()) {

      double distance = Pf7_distance_ptr_->getDistance(x_genome_id, y_genome_id);
      if (std::isnan(distance)) {

        continue;

      } else {

        if (smallest.first > distance) {

          smallest.first = distance;
          smallest.second = x_genome_id;

        }

      }

    } // X

    if (smallest.second == y_genome_id ) {

      ++self_smallest;

    } else {

      ++other_smallest;

    }
    smallest.first = 1000.0;
    smallest.second = "";

  }  // Y

  av_distance = sum_diagonal / static_cast<double>(diagonal_count);
  ExecEnv::log().info("Smallest distance sample,  self smallest: {}, other smallest: {}", self_smallest, other_smallest);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Examine off diagonal
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  sum_diagonal = 0.0;
  diagonal_count = 0;
  nan_count = 0;

  for (size_t y_index = 0; y_index < distance_sample_size; ++y_index) {

    for (size_t x_index = 0; x_index < distance_sample_size; ++x_index) {

      if (x_index == y_index) {

        continue;

      }

      double distance = Pf7_distance_ptr_->getDistance(x_index, y_index);
      if (std::isnan(distance)) {

        ++nan_count;

      } else {

        sum_diagonal += distance;
        ++diagonal_count;

      }

    } // X

  }  // Y

  av_distance = sum_diagonal / static_cast<double>(diagonal_count);
  ExecEnv::log().info("Off Diagonal average distance: {}, Off Diagonal non-nan: {}, Off Diagonal nan: {}", av_distance, diagonal_count, nan_count);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Off diagonal QC Pass filtered
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  sum_diagonal = 0.0;
  diagonal_count = 0;
  nan_count = 0;

  for (auto const& [y_genome_id, y_genome_ptr] : filtered_population_ptr->getMap()) {

    for (auto const& [x_genome_id, x_genome_ptr] : filtered_population_ptr->getMap()) {

      if (y_genome_id == x_genome_id) {

        continue;

      }

      double distance = Pf7_distance_ptr_->getDistance(x_genome_id, y_genome_id);
      if (std::isnan(distance)) {

        ++nan_count;

      } else {

        sum_diagonal += distance;
        ++diagonal_count;

      }

    } // X

  }  // Y

  av_distance = sum_diagonal / static_cast<double>(diagonal_count);
  ExecEnv::log().info("QC Pass; off Diagonal average distance: {}, Off Diagonal non-nan: {}, Off Diagonal nan: {}", av_distance, diagonal_count, nan_count);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Check matrix symmetry.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  size_t sym_equal{0};
  size_t sym_not_equal{0};
  for (size_t y_index = 0; y_index < distance_sample_size; ++y_index) {

    for (size_t x_index = 0; x_index < distance_sample_size; ++x_index) {

      double distance = Pf7_distance_ptr_->getDistance(x_index, y_index);
      double sym_distance = Pf7_distance_ptr_->getDistance(x_index, y_index);
      if (std::isnan(distance) and std::isnan(sym_distance)) {

        ++sym_equal;

      } else if (std::isnan(distance) and not std::isnan(sym_distance)) {

        ++sym_not_equal;

      } else if (not std::isnan(distance) and std::isnan(sym_distance)) {

        ++sym_not_equal;

      } else if (distance == sym_distance) {

        ++sym_equal;

      } else {

        ++sym_not_equal;

      }

    } // X

  }  // Y

  ExecEnv::log().info("Total distance entries: {}, Symmetric: {}, Non Symmetric: {}", (sym_equal + sym_not_equal), sym_equal, sym_not_equal);

}


