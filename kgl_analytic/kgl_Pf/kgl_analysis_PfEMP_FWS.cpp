//
// Created by kellerberrin on 22/05/23.
//

#include "kgl_analysis_PfEMP_FWS.h"
#include "kgl_variant_filter.h"
#include "kgl_analysis_PfEMP_filter.h"
#include <fstream>


namespace kgl = kellerberrin::genome;


void kgl::CalcFWS::calcFwsStatistics(const std::shared_ptr<const PopulationDB>& population) {

// Filter the population variants into frequency bins.

  for (size_t i = 0; i < FWS_FREQUENCY_ARRAY_SIZE; ++i) {

    auto const& [lower_freq, upper_freq] = getFrequency(static_cast<AlleleFrequencyBins>(i));

    P7FrequencyFilter low_freq_filter(FreqInfoField::AF, lower_freq);
    P7FrequencyFilter complement_high_freq_filter(FreqInfoField::AF, upper_freq);
    auto range_filter = AndFilter(low_freq_filter, NotFilter(complement_high_freq_filter));

    std::shared_ptr<const PopulationDB> freq_population = population->viewFilter(range_filter);

    updateGenomeFWSMap(freq_population, i);

  }

}


void kgl::CalcFWS::updateGenomeFWSMap(const std::shared_ptr<const PopulationDB>& freq_population, size_t freq_bin) {

  for (auto const& [genome_id, genome_ptr] : freq_population->getMap()) {

    auto find_iter = genome_fws_map_.find(genome_id);
    if (find_iter == genome_fws_map_.end()) {

      auto [insert_iter, result] = genome_fws_map_.try_emplace(genome_id, FwsFrequencyArray());
      if (not result) {

        ExecEnv::log().error("CalcFWS::updateGenomeFWSMap; Unable to insert genome: {} (duplicate)", genome_id);
        continue;

      }

      find_iter = insert_iter;

    }

    auto& [id, freq_array] = *find_iter;

    auto& freq_record = freq_array[freq_bin];

    updateFreqRecord(genome_ptr, freq_record);

  }

}


void kgl::CalcFWS::updateFreqRecord(const std::shared_ptr<const GenomeDB>& genome_ptr, VariantAnalysisType& freq_record) {

  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [offset, offset_ptr] : contig_ptr->getMap()) {

      HeteroHomoZygous::updateVariantAnalysisType(offset_ptr, freq_record);

    }

  }

}


std::pair<double, double> kgl::CalcFWS::getFrequency(AlleleFrequencyBins bin_type) {

  switch(bin_type) {

    case AlleleFrequencyBins::PERCENT_0_5:
      return {0.0, 0.05 };

    case AlleleFrequencyBins::PERCENT_5_10:
      return {0.05, 0.10 };

    case AlleleFrequencyBins::PERCENT_10_15:
      return {0.10, 0.15 };

    case AlleleFrequencyBins::PERCENT_15_20:
      return {0.15, 0.20 };

    case AlleleFrequencyBins::PERCENT_20_25:
      return {0.20, 0.25 };

    case AlleleFrequencyBins::PERCENT_25_30:
      return {0.25, 0.30 };

    case AlleleFrequencyBins::PERCENT_30_35:
      return {0.30, 0.35 };

    case AlleleFrequencyBins::PERCENT_35_40:
      return {0.35, 0.40 };

    case AlleleFrequencyBins::PERCENT_40_45:
      return {0.40, 0.45 };

    case AlleleFrequencyBins::PERCENT_45_50:
      return {0.45, 0.5 };

  }

  return {0.0 , 0.0}; // Never reached, just to keep the compiler happy

}

void kgl::CalcFWS::writeResults(const std::string& file_name) const {

  std::ofstream analysis_file(file_name);

  if (not analysis_file.good()) {

    ExecEnv::log().error("CalcFWS::writeResults; Unable to open results file: {}", file_name);
    return;

  }

  // Write the header.

  analysis_file << "Genome";

  for (size_t i = 0; i <= FWS_FREQUENCY_ARRAY_SIZE; ++i) {

    analysis_file << CSV_DELIMITER_
                  << "Hom/Het"
                  << CSV_DELIMITER_
                  << "Variant Count"
                  << CSV_DELIMITER_
                  << "Single Variant"
                  << CSV_DELIMITER_
                  << "Homozygous"
                  << CSV_DELIMITER_
                  << "Heterozygous"
                  << CSV_DELIMITER_
                  << "SNP"
                  << CSV_DELIMITER_
                  << "Indel";

  }

  analysis_file << '\n';

  for (auto& [genome_id, freq_array] : genome_fws_map_) {

    analysis_file << genome_id;

    for (size_t i = 0; i <= FWS_FREQUENCY_ARRAY_SIZE; ++i) {

      double hom_het_ratio{0.0};
      size_t het_count = freq_array[i].single_variant_ + freq_array[i].heterozygous_count_;
      if (het_count > 0) {

        hom_het_ratio = static_cast<double>(freq_array[i].homozygous_count_) / static_cast<double>(het_count);

      }

      analysis_file << CSV_DELIMITER_
                    << hom_het_ratio
                    << CSV_DELIMITER_
                    << freq_array[i].total_variants_
                    << CSV_DELIMITER_
                    << freq_array[i].single_variant_
                    << CSV_DELIMITER_
                    << freq_array[i].homozygous_count_
                    << CSV_DELIMITER_
                    << freq_array[i].heterozygous_count_
                    << CSV_DELIMITER_
                    << freq_array[i].snp_count_
                    << CSV_DELIMITER_
                    << freq_array[i].indel_count_;

    }

    analysis_file << '\n';

  }

}
