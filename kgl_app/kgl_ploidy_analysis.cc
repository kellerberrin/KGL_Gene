//
// Created by kellerberrin on 2/03/18.
//


#include "kgl_ploidy_analysis.h"
#include "kgl_lock.h"

#include <fstream>


namespace kgl = kellerberrin::genome;


// This function modifies the data and is thread safe.
bool kgl::PloidyAnalysis::addPloidyRecord(const std::string& genome,
                                          bool homozygous,
                                          bool hq_homozygous,
                                          bool heterozygous,
                                          bool hq_heterozygous,
                                          double ratio) {

  AutoMutex auto_mutex(ploidy_mutex_);

  auto ploidy_record_it = ploidy_data_map_.find(genome);

  if (ploidy_record_it == ploidy_data_map_.end()) {

    std::pair<std::string, PloidyData> inserted_record(genome, PloidyData());
    auto inserted_it = ploidy_data_map_.insert(inserted_record);

    if (not inserted_it.second) {

      ExecEnv::log().error("addPloidyRecord(), could not add ploidy record for genome: {}", genome);
      return false;

    }

    ploidy_record_it = inserted_it.first;

  }

  PloidyData& ploidy_data = ploidy_record_it->second;

  if (homozygous) ++ploidy_data.homozygous_;
  if (hq_homozygous) ++ploidy_data.hq_homozygous_;
  if (heterozygous) ++ploidy_data.heterozygous_;
  if (hq_heterozygous) ++ploidy_data.hq_heterozygous_;
  if (hq_homozygous or hq_heterozygous) {

    ploidy_data.sum_hq_ratio_ += ratio;

    size_t index = ratio * static_cast<size_t>(ratio * PLOIDY_ARRAY_SIZE_);
    if (index >= PLOIDY_ARRAY_SIZE_) {

      index = PLOIDY_ARRAY_SIZE_ - 1;

    }

    PloidyData &ploidy_ratio_data = ratio_array_[index];

    if (homozygous) ++ploidy_ratio_data.homozygous_;
    if (hq_homozygous) ++ploidy_ratio_data.hq_homozygous_;
    if (heterozygous) ++ploidy_ratio_data.heterozygous_;
    if (hq_heterozygous) ++ploidy_ratio_data.hq_heterozygous_;

      ploidy_ratio_data.sum_hq_ratio_ += ratio;

  }

  return true;

}

// Write the data to a results file (not threaded).

bool kgl::PloidyAnalysis::writePloidyResults(const std::string& file_name, char delimiter ) const {

  // open the file.
  std::fstream out_file(file_name, std::fstream::out);
  if (not out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  } else {

    ExecEnv::log().info("Successfully opened file: {} for Ploidy output", file_name);

  }

  out_file << "Genome" << delimiter
           << "Homozygous" << delimiter
           << "HQ_Homozygous" << delimiter
           << "Heterozygous" << delimiter
           << "HQ_Heterozygous" << delimiter
           << "HQ_Ratio" << '\n';

  for( auto ploidy_data : ploidy_data_map_) {

    size_t total_hq = ploidy_data.second.hq_homozygous_ + ploidy_data.second.hq_heterozygous_;
    double hq_ratio;
    if (total_hq != 0) {

      hq_ratio = ploidy_data.second.sum_hq_ratio_ / static_cast<double>(total_hq);

    } else {

      hq_ratio = 0.0;

    }

    out_file << ploidy_data.first << delimiter
             << ploidy_data.second.homozygous_ << delimiter
             << ploidy_data.second.hq_homozygous_ << delimiter
              << ploidy_data.second.heterozygous_ << delimiter
             << ploidy_data.second.hq_heterozygous_ << delimiter
             << hq_ratio << '\n';

  }

  out_file << "Genome" << delimiter
           << "Homozygous" << delimiter
           << "HQ_Homozygous" << delimiter
           << "Heterozygous" << delimiter
           << "HQ_Heterozygous" << delimiter
           << "HQ_Ratio" << '\n';

  size_t count = 0;
  for( auto ploidy_ratio_data : ratio_array_) {

    size_t total_hq = ploidy_ratio_data.hq_homozygous_ + ploidy_ratio_data.hq_heterozygous_;
    double hq_ratio;
    if (total_hq != 0) {

      hq_ratio = ploidy_ratio_data.sum_hq_ratio_ / static_cast<double>(total_hq);

    } else {

      hq_ratio = 0.0;

    }

    count++;
    double bin = static_cast<double>(count) / static_cast<double>(PLOIDY_ARRAY_SIZE_);

    out_file << bin << delimiter
             << ploidy_ratio_data.homozygous_ << delimiter
             << ploidy_ratio_data.hq_homozygous_ << delimiter
             << ploidy_ratio_data.heterozygous_ << delimiter
             << ploidy_ratio_data.hq_heterozygous_ << delimiter
             << hq_ratio << '\n';

  }

  return out_file.good();

}