//
// Created by kellerberrin on 2/03/18.
//


#include "kgl_ploidy_analysis.h"
#include "kgl_lock.h"

#include <fstream>


namespace kgl = kellerberrin::genome;


// This function modifies the data and is thread safe.
bool kgl::PloidyAnalysis::addPloidyRecord(const std::string& genome, bool haploid, bool hq_haploid, bool diploid, bool hq_diploid) {

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

  if (haploid) ++ploidy_data.homozygous_;
  if (hq_haploid) ++ploidy_data.hq_homozygous_;
  if (diploid) ++ploidy_data.heterozygous_;
  if (hq_diploid) ++ploidy_data.hq_heterozygous_;

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
           << "HQ_Heterozygous" << '\n';

  for( auto ploidy_data : ploidy_data_map_) {

    out_file << ploidy_data.first << delimiter
             << ploidy_data.second.homozygous_ << delimiter
             << ploidy_data.second.hq_homozygous_ << delimiter
             << ploidy_data.second.heterozygous_ << delimiter
             << ploidy_data.second.hq_heterozygous_ << '\n';

  }

  return out_file.good();

}