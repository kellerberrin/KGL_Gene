//
// Created by kellerberrin on 27/04/23.
//


#include "kgl_variant_filter_db_variant.h"
#include "kgl_variant_filter_Pf7.h"
#include "kga_analysis_lib_PfFilter.h"


namespace kga = kellerberrin::genome::analysis;
namespace kgl = kellerberrin::genome;


// Quality filter the variants using read depth, VQSLOD and other statistics
std::shared_ptr<kgl::PopulationDB> kga::FilterPf7::qualityFilter(const std::shared_ptr<const PopulationDB>& unfiltered_population_ptr) const {


  size_t unfiltered_count = unfiltered_population_ptr->variantCount();
  ExecEnv::log().info("Unfiltered Population: {}, Genome count: {}, Variant Count: {}",
                      unfiltered_population_ptr->populationId(),
                      unfiltered_population_ptr->getMap().size(),
                      unfiltered_count);


  auto filtered_qc_population_ptr = unfiltered_population_ptr;
  size_t qc_count{0};
  if (filter_qc_active_) {
    // Shallow filter only.
    filtered_qc_population_ptr = Pf7_sample_ptr_->filterPassQCGenomes(unfiltered_population_ptr);

    qc_count = filtered_qc_population_ptr->variantCount();
    double qc_filtered = 100.0 * (static_cast<double>(unfiltered_count - qc_count) / static_cast<double>(unfiltered_count));
    ExecEnv::log().info("Filter P7 QC Pass Genome count: {}, Variant Count: {}, Sample Data Count: {}, filtered: {:.2f}%",
                        filtered_qc_population_ptr->getMap().size(),
                        qc_count,
                        Pf7_sample_ptr_->getMap().size(),
                        qc_filtered);

  } else {

    qc_count = filtered_qc_population_ptr->variantCount();

  }

  auto monoclonal_population_ptr = filtered_qc_population_ptr;
  if (filter_fws_active_) {

    // Shallow filter only.
    monoclonal_population_ptr = Pf7_fws_ptr_->viewFilterFWS(FwsFilterType::GREATER_EQUAL, fws_monoclonal_threshold_, filtered_qc_population_ptr);

    size_t monoclonal_count = monoclonal_population_ptr->variantCount();
    double mono_filtered = 100.0 * (static_cast<double>(qc_count - monoclonal_count) / static_cast<double>(qc_count));
    ExecEnv::log().info("Filter MonoClonal FWS: {}, Genome Count: {},Variant Count: {}, filtered: {:.2f}%",
                        fws_monoclonal_threshold_,
                        monoclonal_population_ptr->getMap().size(),
                        monoclonal_count,
                        mono_filtered);
  }

  // Call bespoke variant filter.
  P7VariantFilter info_field_filter;
  P7VariantFilter::initializeStats();
  auto filtered_population_ptr = monoclonal_population_ptr->viewFilter(info_field_filter);
  P7VariantFilter::printStats();

  // AF Frequency Filter.
  if (af_filter_active_) {

    P7FrequencyFilter af_frequency_filter(FreqInfoField::AF, variant_frequency_cutoff_);
    P7FrequencyFilter::initializeStats();
    filtered_population_ptr = filtered_population_ptr->viewFilter(af_frequency_filter);
    P7FrequencyFilter::printStats(variant_frequency_cutoff_);

  }

  // MLEAF Frequency Filter (redundant).
  if (mleaf_filter_active_) {

    P7FrequencyFilter mleaf_frequency_filter(FreqInfoField::MLEAF, variant_frequency_cutoff_);
    P7FrequencyFilter::initializeStats();
    filtered_population_ptr = filtered_population_ptr->viewFilter(mleaf_frequency_filter);
    P7FrequencyFilter::printStats(variant_frequency_cutoff_);

  }

  // Filter for snp only.
  if (snp_filter_active_) {

    filtered_population_ptr = filtered_population_ptr->viewFilter(SNPFilter());

  }

  if (coding_filter_active_) {

// Recode the coding variants filter using gene, transcript functionality.
//    filtered_population_ptr = filtered_population_ptr->viewFilter(FilterAllCodingVariants(genome_3D7_ptr_));
    ExecEnv::log().info("** Coding Filter Not Active** Coding Population Final Filtered Size Genome count: {}, Variant Count: {}",
                        filtered_population_ptr->getMap().size(),
                        filtered_population_ptr->variantCount());


  }

  // We need to do a deep copy of the filtered population here since the pass QC and FWS P7 filters only do a shallow copy.
  // And when the resultant population pointers go out of scope they will take the shared population structure with them.
  auto deepcopy_population_ptr = filtered_population_ptr->deepCopy();

  // Filtered population should contain all contigs for all genomes.
  deepcopy_population_ptr->squareContigs();

  ExecEnv::log().info("Population Final Filtered Size Genome count: {}, Variant Count: {}",
                      deepcopy_population_ptr->getMap().size(),
                      deepcopy_population_ptr->variantCount());

  return deepcopy_population_ptr;

}



// Population must be PF7
// Important - this code above only filters a shallow copy of the population.
std::shared_ptr<kgl::PopulationDB> kgl::Pf7SampleResource::filterPassQCGenomes(const std::shared_ptr<const PopulationDB>& Pf7_unfiltered_ptr) const {


  auto filtered_ptr = std::make_shared<kgl::PopulationDB>(Pf7_unfiltered_ptr->populationId() + "_QC_Pass",
                                                          Pf7_unfiltered_ptr->dataSource());

  for (auto const& [genome_id, genome_ptr] : Pf7_unfiltered_ptr->getMap()) {

    if (getMap().contains(genome_id)) {

      auto sample_iter = getMap().find(genome_id);
      auto const& [id, sample_data] = *sample_iter;
      if (sample_data.pass()) {

        if (not filtered_ptr->addGenome(genome_ptr)) {

          ExecEnv::log().warn("Pf7SampleResource::filterPassQCGenomes; Unable to add filtered genome: {} to filtered  population", genome_id);

        }

      } // Pass

    } else {

      ExecEnv::log().warn("PPf7SampleResource::filterPassQCGenomes; Genome: {} not found in sample data", genome_id);

    }

  } // For genomes.

  return filtered_ptr;

}


// Important - this code above only filters a shallow copy of the population.
std::shared_ptr<kgl::PopulationDB> kga::FilterPf3k::filterCOI(const std::shared_ptr<const PopulationDB>& population_ptr) const {

  auto filtered_ptr = std::make_shared<kgl::PopulationDB>(population_ptr->populationId(), population_ptr->dataSource());

  size_t no_COI_record{0};
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    auto coi_opt = Pf3k_COI_ptr_->genomeCOI(genome_id);
    if (coi_opt) {

      auto coi_value = coi_opt.value();
      if (coi_value == 1) {

        if (not filtered_ptr->addGenome(genome_ptr)) {

          ExecEnv::log().warn("FilterPf3k::filterCOI; Unable to add filtered genome: {} to filtered  population", genome_id);

        }

      } // Pass

    } else {

      ++no_COI_record;

    }

  } // For genomes.

  ExecEnv::log().info("FilterPf3k::filterCOI; Unfiltered genomes: {}, COI = 1 genomes: {}, No COI record genomes: {}",
                      population_ptr->getMap().size(),
                      filtered_ptr->getMap().size(),
                      no_COI_record);

  return filtered_ptr;

}

