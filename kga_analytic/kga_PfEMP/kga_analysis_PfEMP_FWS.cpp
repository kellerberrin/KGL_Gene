//
// Created by kellerberrin on 22/05/23.
//

#include "kga_analysis_PfEMP_FWS.h"
#include "kgl_variant_filter_db_variant.h"
#include "kgl_variant_filter_Pf7.h"
#include "kgl_variant_db_variant.h"
#include <fstream>


namespace kgl = kellerberrin::genome;


void kgl::CalcFWS::calcFwsStatistics(const std::shared_ptr<const PopulationDB>& population) {

  ExecEnv::log().info("CalcFWS::calcFwsStatistics; total population variants: {}", population->variantCount());

// Calc variant statistics.
  updateVariantFWSMap(population);

// Filter the population variants into frequency bins.
  for (size_t i = 0; i < FWS_FREQUENCY_ARRAY_SIZE; ++i) {

    auto const& [lower_freq, upper_freq] = getFrequency(static_cast<AlleleFrequencyBins>(i));

    P7FrequencyFilter low_freq_filter(FreqInfoField::AF, lower_freq);
    P7FrequencyFilter complement_high_freq_filter(FreqInfoField::AF, upper_freq);
    auto range_filter = AndFilter(low_freq_filter, NotFilter(complement_high_freq_filter));

    std::shared_ptr<const PopulationDB> freq_population = population->viewFilter(range_filter);

    ExecEnv::log().info("CalcFWS::calcFwsStatistics; lower: {}, upper: {}, population variants: {}", lower_freq, upper_freq, freq_population->variantCount());

    updateGenomeFWSMap(freq_population, i);

  }

}

void kgl::CalcFWS::updateVariantFWSMap(const std::shared_ptr<const PopulationDB>& population) {

  VariantDBVariant variant_db_variant(population);

  for (auto const& [hgvs, variant_record] : variant_db_variant.variantMap()) {

    auto find_iter = variant_fws_map_.find(hgvs);
    if (find_iter == variant_fws_map_.end()) {

      auto [insert_iter, result] = variant_fws_map_.try_emplace(hgvs, AlleleSummmary());
      if (not result) {

        ExecEnv::log().error("CalcFWS::updateVariantFWSMap; Unable to insert Variant: {} (duplicate)", hgvs);
        continue;

      }

      find_iter = insert_iter;

    }

    auto const& [variant_ptr, variant_index] = variant_record;
    auto variant_summary = variant_db_variant.summaryByVariant(variant_ptr);
    auto& [map_hgvs, map_summary] = *find_iter;

    map_summary += variant_summary;

  }

}

void kgl::CalcFWS::updateGenomeFWSMap(const std::shared_ptr<const PopulationDB>& freq_population, size_t freq_bin) {

  VariantDBVariant variant_db_variant(freq_population);

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

    auto genome_summary = variant_db_variant.summaryByGenome(genome_id);
    freq_record += genome_summary;

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

    case AlleleFrequencyBins::PERCENT_50_100:
      return {0.5, 1.0 };

  }

  return {0.0 , 0.0}; // Never reached, just to keep the compiler happy

}

void kgl::CalcFWS::writeGenomeResults(const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr, const std::string& file_name) const {

  std::ofstream analysis_file(file_name);

  if (not analysis_file.good()) {

    ExecEnv::log().error("CalcFWS::writeGenomeResults; Unable to open results file: {}", file_name);
    return;

  }

  // Write the header.

  analysis_file << "Genome"
                << CSV_DELIMITER_
                << "FWS";

  for (size_t i = 0; i < FWS_FREQUENCY_ARRAY_SIZE; ++i) {

    analysis_file << CSV_DELIMITER_
                  << "LowerFreq"
                  << CSV_DELIMITER_
                  << "UpperFreq"
                  << CSV_DELIMITER_
                  << "Hom/Het"
                  << CSV_DELIMITER_
                  << "Minor Hom/Het"
                  << CSV_DELIMITER_
                  << "Variant Count"
                  << CSV_DELIMITER_
                  << "Hom Ref (A;A)"
                  << CSV_DELIMITER_
                  << "Het Ref Minor (A;a)"
                  << CSV_DELIMITER_
                  << "Hom Minor (a;a)";

  }

  analysis_file << '\n';

  for (auto& [genome_id, freq_array] : genome_fws_map_) {

    analysis_file << genome_id
                  << CSV_DELIMITER_
                  << Pf7_fws_ptr->getFWS(genome_id);

    for (size_t i = 0; i < FWS_FREQUENCY_ARRAY_SIZE; ++i) {

      double hom_het_ratio{0.0};
      size_t het_count = freq_array[i].minorHeterozygous_;
      size_t hom_count = freq_array[i].minorHomozygous_ + freq_array[i].referenceHomozygous_;
      if (het_count > 0) {

        hom_het_ratio = static_cast<double>(hom_count) / static_cast<double>(het_count);

      }

      double minor_hom_het_ratio{0.0};
      size_t minor_het_count = freq_array[i].minorHeterozygous_;
      size_t minor_hom_count = freq_array[i].minorHomozygous_;
      if (minor_het_count > 0) {

        minor_hom_het_ratio = static_cast<double>(minor_hom_count) / static_cast<double>(minor_het_count);

      }

      size_t total_variants = freq_array[i].minorHeterozygous_ + freq_array[i].minorHomozygous_;

      auto const& [lower_range, upper_range] = getFrequency(static_cast<AlleleFrequencyBins>(i));

      analysis_file << CSV_DELIMITER_
                    << lower_range
                    << CSV_DELIMITER_
                    << upper_range
                    << CSV_DELIMITER_
                    << hom_het_ratio
                    << CSV_DELIMITER_
                    << minor_hom_het_ratio
                    << CSV_DELIMITER_
                    << total_variants
                    << CSV_DELIMITER_
                    << freq_array[i].referenceHomozygous_
                    << CSV_DELIMITER_
                    << freq_array[i].minorHeterozygous_
                    << CSV_DELIMITER_
                    << freq_array[i].minorHomozygous_;

    }

    analysis_file << '\n';

  }

}


void kgl::CalcFWS::writeVariantResults(const std::string& file_name) const {

  std::ofstream analysis_file(file_name);

  if (not analysis_file.good()) {

    ExecEnv::log().error("CalcFWS::writeVariantResults; Unable to open results file: {}", file_name);
    return;

  }

  // Write the header.
   analysis_file << "Variant"
                << CSV_DELIMITER_
                << "Hom/Het"
                << CSV_DELIMITER_
                << "Minor Hom/Het"
                << CSV_DELIMITER_
                << "Genome Count"
                << CSV_DELIMITER_
                << "Hom Ref (A;A)"
                << CSV_DELIMITER_
                << "Het Ref Minor (A;a)"
                << CSV_DELIMITER_
                << "Hom Minor (a;a)"
                << '\n';

  for (auto& [hgvs, summary] : variant_fws_map_) {

    double hom_het_ratio{0.0};
    size_t het_count = summary.minorHeterozygous_;
    size_t hom_count = summary.minorHomozygous_ + summary.referenceHomozygous_;
    if (het_count > 0) {

      hom_het_ratio = static_cast<double>(hom_count) / static_cast<double>(het_count);

    }

    double minor_hom_het_ratio{0.0};
    size_t minor_het_count = summary.minorHeterozygous_;
    size_t minor_hom_count = summary.minorHomozygous_;
    if (minor_het_count > 0) {

      minor_hom_het_ratio = static_cast<double>(minor_hom_count) / static_cast<double>(minor_het_count);

    }

    size_t total_genomes = summary.minorHeterozygous_ + summary.minorHomozygous_;

    analysis_file << hgvs
                  << CSV_DELIMITER_
                  << hom_het_ratio
                  << CSV_DELIMITER_
                  << minor_hom_het_ratio
                  << CSV_DELIMITER_
                  << total_genomes
                  << CSV_DELIMITER_
                  << summary.referenceHomozygous_
                  << CSV_DELIMITER_
                  << summary.minorHeterozygous_
                  << CSV_DELIMITER_
                  << summary.minorHomozygous_
                  << '\n';

  }

}
