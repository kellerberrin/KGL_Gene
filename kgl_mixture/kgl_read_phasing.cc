//
// Created by kellerberrin on 27/08/18.
//

#include "kgl_read_phasing.h"
#include "kel_exec_env.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>


namespace kgl = kellerberrin::genome;
namespace bt = boost;



bool kgl::GenomeMixtureStatistics::getMixtureStatistics(const GenomeId_t& genome,
                                                        MixtureStatistics& mixture_statistics) const {


  auto result = getMap().find(genome);

  if (result == getMap().end()) {

    mixture_statistics.first = 0;
    for (auto& strain : mixture_statistics.second) {

      strain = 0.0;

    }

    return false;

  }

  mixture_statistics = result->second;
  return true;

}



bool kgl::GenomeMixtureStatistics::readMixtureStatistics(const std::string& file_name, const std::string& field_delimiters) {

  std::ifstream readfile(file_name);

  if (not readfile.good()) {

    ExecEnv::log().error("GenomeMixtureStatistics::readMixtureStatistics; could not open mixture file: {}", file_name);
    return false;

  }

  bt::char_separator<char> item_key_sep(field_delimiters.c_str());
  size_t line_count = 0;

  while (true) {

    std::string record_str;
    if (std::getline(readfile, record_str).eof()) break;

    if ((record_str)[0] == COMMENT_CHAR_) continue;   // Ignore comment lines.

    ++line_count;
    std::vector<std::string> field_vector;
    bt::tokenizer<bt::char_separator<char>> tokenize_item(record_str, item_key_sep);

    for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

      field_vector.push_back(*iter_item);

    }

    if (field_vector.size() < MIN_FIELD_COUNT_) {

      ExecEnv::log().error("GenomeMixtureStatistics::readMixtureStatistics; Required minimum fields : {}, actual fields: {}, in line: {} text: {}",
                           MIN_FIELD_COUNT_, field_vector.size(), line_count, record_str);
      continue;

    }

    try {

      double p1 = std::stof(field_vector[1]);
      double max = p1;
      double p2 = std::stof(field_vector[2]);
      max = max < p2 ? p2 : max;
      double p3 = std::stof(field_vector[3]);
      max = max < p3 ? p3 : max;
      double p4 = std::stof(field_vector[4]);
      max = max < p4 ? p4 : max;
      double p5 = std::stof(field_vector[5]);
      max = max < p5 ? p5 : max;
      size_t strain = std::stoll(field_vector[6]);

      if (std::fabs(1.0-(p1+p2+p3+p4+p5)) > TOLERANCE_) {

        ExecEnv::log().warn("GenomeMixtureStatistics::readMixtureStatistics; Strain proportions: {} do not sum to 1.0 in line: {} text: {}",
                             (p1+p2+p3+p4+p5), line_count, record_str);
        continue;

      }

      if (strain == 1 and max < SINGLE_STRAIN_) {

        ExecEnv::log().warn("GenomeMixtureStatistics::readMixtureStatistics; Single Strain proportion: {} less than : {} in line: {} text: {}",
                             max, SINGLE_STRAIN_, line_count, record_str);
        continue;

      }

      MixtureStatistics mixture_statistics;
      mixture_statistics.first = strain;
      mixture_statistics.second[0] = p1;
      mixture_statistics.second[1] = p2;
      mixture_statistics.second[2] = p3;
      mixture_statistics.second[3] = p4;
      mixture_statistics.second[4] = p5;

      std::pair<GenomeId_t, MixtureStatistics> insert_record(field_vector[0], mixture_statistics);
      auto result = mixture_statistics_map_.insert(insert_record);

      if (not result.second) {

        ExecEnv::log().error("GenomeMixtureStatistics::readMixtureStatistics; Unable to add mixture record for genome: {} (duplicate) for line: {} text: {}",
                             field_vector[0], line_count, record_str);

      }

    } catch(...) {

      ExecEnv::log().error("GenomeMixtureStatistics::readMixtureStatistics; Unexpected numeric field format in line: {} text: {}", line_count, record_str);

    }

  }

  readfile.close();
  ExecEnv::log().info("GenomeMixtureStatistics::readMixtureStatistics; Mixture file: {}, read: {} genomes", file_name, line_count);
  return true;

}

