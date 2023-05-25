//
// Created by kellerberrin on 24/05/23.
//

#include "kgl_variant_filter_Pf7.h"
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

