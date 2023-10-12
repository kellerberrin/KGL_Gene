//
// Created by kellerberrin on 27/04/23.
//


#include "kgl_variant_filter_db_variant.h"
#include "kgl_variant_filter_Pf7.h"
#include "kgl_analysis_PfEMP.h"


namespace kgl = kellerberrin::genome;


// Quality filter the variants using read depth, VQSLOD and other statistics
std::shared_ptr<kgl::PopulationDB> kgl::PfEMPAnalysis::qualityFilter(const std::shared_ptr<const PopulationDB>& unfiltered_population_ptr) {


  size_t unfiltered_count = unfiltered_population_ptr->variantCount();
  ExecEnv::log().info("Unfiltered Population: {}, Genome count: {}, Variant Count: {}",
                      unfiltered_population_ptr->populationId(),
                      unfiltered_population_ptr->getMap().size(),
                      unfiltered_count);


  auto filtered_qc_population_ptr = unfiltered_population_ptr;
  size_t qc_count{0};
  if constexpr (FILTER_QC_ACTIVE_) {
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
  if constexpr (FILTER_FWS_ACTIVE_) {

    // Shallow filter only.
    monoclonal_population_ptr = Pf7_fws_ptr_->viewFilterFWS(FwsFilterType::GREATER_EQUAL, FWS_MONOCLONAL_THRESHOLD, filtered_qc_population_ptr);

    size_t monoclonal_count = monoclonal_population_ptr->variantCount();
    double mono_filtered = 100.0 * (static_cast<double>(qc_count - monoclonal_count) / static_cast<double>(qc_count));
    ExecEnv::log().info("Filter MonoClonal FWS: {}, Genome Count: {},Variant Count: {}, filtered: {:.2f}%",
                        FWS_MONOCLONAL_THRESHOLD,
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
  if constexpr (AF_FILTER_ACTIVE_) {

    P7FrequencyFilter af_frequency_filter(FreqInfoField::AF, VARIANT_FREQUENCY_CUTOFF_);
    P7FrequencyFilter::initializeStats();
    filtered_population_ptr = filtered_population_ptr->viewFilter(af_frequency_filter);
    P7FrequencyFilter::printStats(VARIANT_FREQUENCY_CUTOFF_);

  }

  // MLEAF Frequency Filter (redundant).
  if constexpr (MLEAF_FILTER_ACTIVE_) {

    P7FrequencyFilter mleaf_frequency_filter(FreqInfoField::MLEAF, VARIANT_FREQUENCY_CUTOFF_);
    P7FrequencyFilter::initializeStats();
    filtered_population_ptr = filtered_population_ptr->viewFilter(mleaf_frequency_filter);
    P7FrequencyFilter::printStats(VARIANT_FREQUENCY_CUTOFF_);

  }

  // Filter for snp only.
  if constexpr (SNP_FILTER_ACTIVE_) {

    filtered_population_ptr = filtered_population_ptr->viewFilter(SNPFilter());

  }

  if constexpr(CODING_FILTER_ACTIVE_) {

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

