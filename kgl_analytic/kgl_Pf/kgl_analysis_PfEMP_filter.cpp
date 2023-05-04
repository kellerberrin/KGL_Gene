//
// Created by kellerberrin on 27/04/23.
//

#include "kgl_analysis_PfEMP_filter.h"
#include "kgl_analysis_PfEMP.h"
#include "kgl_variant_filter_info.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Frequency filter using AF or MLEAF
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::P7FrequencyFilter::applyFilter(const Variant &variant) const {

  size_t alt_count = variant.evidence().altVariantCount();
  size_t alt_index = variant.evidence().altVariantIndex();

  ++freq_filter_stats_.unfiltered_;

  std::string info_field = info_field_ == FreqInfoField::AF ? AF_FIELD_ : MLEAF_FIELD_;
  auto info_opt = InfoEvidenceAnalysis::getTypedInfoData<std::vector<double>>(variant, info_field);
  if (info_opt) {

    std::vector<double> info_vector = std::move(info_opt.value());
    if (info_vector.size() != alt_count) {

      ExecEnv::log().error("P7FrequencyFilter::applyFilter; AF vector size: {}, not equal alt variant count: {}, Info field: {}",
                           info_vector.size(), alt_count, info_field);

      return false;

    }
    if (info_vector.size() <= alt_index) {

      ExecEnv::log().error("P7FrequencyFilter::applyFilter; alt variant index: {} out of range for vector size:{}, Info field: {}",
                           alt_index, info_vector.size(), info_field);
      return false;

    }

    if (info_vector[alt_index] >= freq_cutoff_) {

      ++freq_filter_stats_.accepted_;
      return true;

    } else {

      ++freq_filter_stats_.rejected_;
      return false;

    }

  }

  ++freq_filter_stats_.missing_;

  return true;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simple object to collect filter statistics for each filter field.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::FilterFieldInfo::printStats(double filter_level) const {

  ExecEnv::log().info("Filter: {}, Level: {}, Variants: {}, Accepted: {:.2f}%, Rejected: {:.2f}%, Rejection Rate: {:.2f}%, Missing Info Field: {:.2f}%",
                      filter_name_,
                      filter_level,
                      unfiltered_.load(),
                      acceptedPercent(),
                      rejectedPercent(),
                      rejectionRate(),
                      missingPercent());

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Bespoke variant quality filter for the P7 Pf database.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Before filtering.
void kgl::P7VariantFilter::initializeStats() {

  vqslod_stats_.initialize();
  qd_stats_.initialize();
  mq_stats_.initialize();
  sor_stats_.initialize();
  mqrs_stats_.initialize();
  rprs_stats_.initialize();
  if constexpr (READDEPTH_ACTIVE_) depth_stats_.initialize();
  unfiltered_variants_ = 0;
  accepted_variants_ = 0;

}

// After filtering.
void kgl::P7VariantFilter::printStats() {

  vqslod_stats_.printStats(VQSLOD_LEVEL_);
  qd_stats_.printStats(QD_LEVEL_);
  mq_stats_.printStats(MQ_LEVEL_),
  sor_stats_.printStats(SOR_LEVEL_);
  mqrs_stats_.printStats(MQRANKSUM_LEVEL_);
  rprs_stats_.printStats(READPOSRANKSUM_LEVEL_);
  if constexpr (READDEPTH_ACTIVE_)  depth_stats_.printStats(MINIMUM_READDEPTH_);

  double accepted = unfiltered_variants_ > 0 ? (100.0 * (static_cast<double>(accepted_variants_) / static_cast<double>(unfiltered_variants_))) : 0.0;
  double rejected = unfiltered_variants_ > 0 ? (100.0 - accepted) : 0.0;
  ExecEnv::log().info("Total UnFiltered Variants: {}, Total Filtered_variants: {}, Accepted: {:.2f}%, Rejected: {:.2f}%",
                      unfiltered_variants_.load(), accepted_variants_.load(), accepted, rejected );

  ExecEnv::log().info("The 'Pf7_Pf3D7_MITO' contig 'MQ' filter accept threshold: {}", Pf7_Pf3D7_MITO_MQ_LEVEL_);

}


// Bespoke general filter for the P7 database.
bool kgl::P7VariantFilter::applyFilter(const Variant & variant) const {

  ++unfiltered_variants_;

// if VQSLOD is present then only use this filter
////////////////////////////////////////////////////////////////////////////////////////////////////
  ++vqslod_stats_.unfiltered_;
  auto info_opt = InfoEvidenceAnalysis::getTypedInfoData<double>( variant, VQSLOD_FIELD_);
  if (info_opt) {

    if (info_opt.value() >= VQSLOD_LEVEL_) {

      ++vqslod_stats_.accepted_;
      ++accepted_variants_;
      return true;

    } else {

      ++vqslod_stats_.rejected_;
      return false;

    }

  }
  ++vqslod_stats_.missing_;

  // If the variant is missing the VQSLOD field, then use all the following filters.
  // QD
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ++qd_stats_.unfiltered_;
  info_opt = InfoEvidenceAnalysis::getTypedInfoData<double>( variant, QD_FIELD_);
  if (info_opt) {

    if (info_opt.value() >= QD_LEVEL_) {

      ++qd_stats_.accepted_;

    } else {

      ++qd_stats_.rejected_;
      return false;

    }

  } else {

    ++qd_stats_.missing_;

  }

  // MQ
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ++mq_stats_.unfiltered_;
  info_opt = InfoEvidenceAnalysis::getTypedInfoData<double>( variant, MQ_FIELD_);
  if (info_opt) {

    // Variants in the 'Pf7_Pf3D7_MITO' contig have a lower MQ threshold.
    double mq_level = variant.contigId() == Pf7_Pf3D7_MITO ? Pf7_Pf3D7_MITO_MQ_LEVEL_ : MQ_LEVEL_;
    if (info_opt.value() >= mq_level) {

      ++mq_stats_.accepted_;

    } else {

      ++mq_stats_.rejected_;
      return false;

    }

  } else {

    ++mq_stats_.missing_;

  }

  // SOR
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ++sor_stats_.unfiltered_;
  info_opt = InfoEvidenceAnalysis::getTypedInfoData<double>( variant, SOR_FIELD_);
  if (info_opt) {

    if (info_opt.value() <= SOR_LEVEL_) {

      ++sor_stats_.accepted_;

    } else {

      ++sor_stats_.rejected_;
      return false;

    }

  } else {

    ++sor_stats_.missing_;

  }

  // MQRANKSUM
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ++mqrs_stats_.unfiltered_;
  info_opt = InfoEvidenceAnalysis::getTypedInfoData<double>( variant, MQRANKSUM_FIELD_);
  if (info_opt) {

    if (info_opt.value() >= MQRANKSUM_LEVEL_) {

      ++mqrs_stats_.accepted_;

    } else {

      ++mqrs_stats_.rejected_;
      return false;

    }

  } else {

    ++mqrs_stats_.missing_;

  }

  // READPOSRANKSUM
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ++rprs_stats_.unfiltered_;
  info_opt = InfoEvidenceAnalysis::getTypedInfoData<double>( variant, READPOSRANKSUM_FIELD_);
  if (info_opt) {

    if (info_opt.value() >= READPOSRANKSUM_LEVEL_) {

      ++rprs_stats_.accepted_;

    } else {

      ++rprs_stats_.rejected_;
      return false;

    }

  } else {

    ++rprs_stats_.missing_;

  }

  // READDEPTH
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  if constexpr (READDEPTH_ACTIVE_) {

    ++depth_stats_.unfiltered_;
    if (variant.evidence().formatData()) {

      auto const &format_data = *(variant.evidence().formatData().value());
      if (format_data.DPCount() >= MINIMUM_READDEPTH_) {

        ++depth_stats_.accepted_;

      } else {

        ++depth_stats_.rejected_;
        return false;

      }

    } else {

      ++depth_stats_.missing_;

    }

  } // if constexpr (READDEPTH_ACTIVE_)

  ++accepted_variants_;

  return true;

}


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
    monoclonal_population_ptr = Pf7_fws_ptr_->filterFWS(FwsFilterType::GREATER_EQUAL, MONOCLONAL_FWS_THRESHOLD, filtered_qc_population_ptr);

    size_t monoclonal_count = monoclonal_population_ptr->variantCount();
    double mono_filtered = 100.0 * (static_cast<double>(qc_count - monoclonal_count) / static_cast<double>(qc_count));
    ExecEnv::log().info("Filter MonoClonal FWS: {}, Genome Count: {},Variant Count: {}, filtered: {:.2f}%",
                        MONOCLONAL_FWS_THRESHOLD,
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
  P7FrequencyFilter af_frequency_filter(FreqInfoField::AF, VARIANT_FREQUENCY_CUTOFF_);
  P7FrequencyFilter::initializeStats();
  filtered_population_ptr  = filtered_population_ptr->viewFilter(af_frequency_filter);
  P7FrequencyFilter::printStats(VARIANT_FREQUENCY_CUTOFF_);

  // MLEAF Frequency Filter (redundant).
  if constexpr (MLEAF_FILTER_ACTIVE_) {

    P7FrequencyFilter mleaf_frequency_filter(FreqInfoField::MLEAF, VARIANT_FREQUENCY_CUTOFF_);
    P7FrequencyFilter::initializeStats();
    filtered_population_ptr = filtered_population_ptr->viewFilter(mleaf_frequency_filter);
    P7FrequencyFilter::printStats(VARIANT_FREQUENCY_CUTOFF_);

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

