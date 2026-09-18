//
// Created by kellerberrin on 24/05/23.
//

#include "kgl_variant_filter_Pf7.h"
#include "kgl_variant_filter_info.h"


namespace kgl = kellerberrin::genome;


/// Frequency filter using AF or MLEAF.
bool kgl::P7FrequencyFilter::applyFilter(const Variant &variant) const {

  size_t alt_count = variant.evidence().altVariantCount();
  size_t alt_index = variant.evidence().altVariantIndex();

  freq_filter_stats_.countUnfiltered();

  auto info_opt = InfoEvidenceAnalysis::getTypedInfoData<std::vector<double>>(variant, info_field_);
  if (info_opt) {

    std::vector<double> info_vector = std::move(info_opt.value());
    if (info_vector.size() != alt_count) {

      ExecEnv::log().error("P7FrequencyFilter::applyFilter; AF vector size: {}, not equal alt variant count: {}, Info field: {}",
                           info_vector.size(), alt_count, info_field_);

      return false;

    }
    if (info_vector.size() <= alt_index) {

      ExecEnv::log().error("P7FrequencyFilter::applyFilter; alt variant index: {} out of range for vector size:{}, Info field: {}",
                           alt_index, info_vector.size(), info_field_);
      return false;

    }

    if (info_vector[alt_index] >= freq_cutoff_) {

      freq_filter_stats_.countAccepted();
      return true;

    }

    freq_filter_stats_.countRejected();
    return false;

  }

  freq_filter_stats_.countMissing();

  return true;

}


/// Simple object to collect filter statistics for each filter field.
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


/// Checks a double info field against a threshold; pass_if_greater_equal selects >= vs <=.
kgl::P7VariantFilter::FieldOutcome kgl::P7VariantFilter::checkInfoField(const Variant& variant, const char* field, double level, bool pass_if_greater_equal) {

  auto info_opt = InfoEvidenceAnalysis::getTypedInfoData<double>(variant, field);
  if (not info_opt) {

    return FieldOutcome::MISSING;

  }

  bool pass = pass_if_greater_equal ? info_opt.value() >= level : info_opt.value() <= level;
  return pass ? FieldOutcome::PASS : FieldOutcome::FAIL;

}


/// Bespoke variant quality filter for the P7 Pf database.
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
  mq_stats_.printStats(MQ_LEVEL_);
  sor_stats_.printStats(SOR_LEVEL_);
  mqrs_stats_.printStats(MQRANKSUM_LEVEL_);
  rprs_stats_.printStats(READPOSRANKSUM_LEVEL_);
  if constexpr (READDEPTH_ACTIVE_)  depth_stats_.printStats(MINIMUM_READDEPTH_);

  double accepted = unfiltered_variants_ > 0 ? (100.0 * (static_cast<double>(accepted_variants_) / static_cast<double>(unfiltered_variants_))) : 0.0;
  double rejected = unfiltered_variants_ > 0 ? (100.0 - accepted) : 0.0;
  ExecEnv::log().info("Total UnFiltered Variants: {}, Total Filtered_variants: {}, Accepted: {:.2f}%, Rejected: {:.2f}%",
                      unfiltered_variants_.load(), accepted_variants_.load(), accepted, rejected );

  ExecEnv::log().info("The 'Pf7_Pf3D7_MITO' contig_ref_ptr 'MQ' filter accept threshold: {}", Pf7_Pf3D7_MITO_MQ_LEVEL_);

}


// Bespoke general filter for the P7 database.
bool kgl::P7VariantFilter::applyFilter(const Variant & variant) const {

  ++unfiltered_variants_;

  // if VQSLOD is present then only use this filter
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  vqslod_stats_.countUnfiltered();
  switch (checkInfoField(variant, VQSLOD_FIELD_, VQSLOD_LEVEL_, true)) {

    case FieldOutcome::PASS:
      vqslod_stats_.countAccepted();
      ++accepted_variants_;
      return true;

    case FieldOutcome::FAIL:
      vqslod_stats_.countRejected();
      return false;

    case FieldOutcome::MISSING:
      vqslod_stats_.countMissing();
      break;

  }

  // If the variant is missing the VQSLOD field, then use all the following filters.
  // QD
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  qd_stats_.countUnfiltered();
  switch (checkInfoField(variant, QD_FIELD_, QD_LEVEL_, true)) {

    case FieldOutcome::MISSING: qd_stats_.countMissing();  break;
    case FieldOutcome::PASS:    qd_stats_.countAccepted(); break;
    case FieldOutcome::FAIL:    qd_stats_.countRejected(); return false;

  }

  // MQ
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  mq_stats_.countUnfiltered();
  // Variants in the 'Pf7_Pf3D7_MITO' contig_ref_ptr have a lower MQ threshold.
  double mq_level = variant.contigId() == Pf7_Pf3D7_MITO ? Pf7_Pf3D7_MITO_MQ_LEVEL_ : MQ_LEVEL_;
  switch (checkInfoField(variant, MQ_FIELD_, mq_level, true)) {

    case FieldOutcome::MISSING: mq_stats_.countMissing();  break;
    case FieldOutcome::PASS:    mq_stats_.countAccepted(); break;
    case FieldOutcome::FAIL:    mq_stats_.countRejected(); return false;

  }

  // SOR
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  sor_stats_.countUnfiltered();
  switch (checkInfoField(variant, SOR_FIELD_, SOR_LEVEL_, false)) {

    case FieldOutcome::MISSING: sor_stats_.countMissing();  break;
    case FieldOutcome::PASS:    sor_stats_.countAccepted(); break;
    case FieldOutcome::FAIL:    sor_stats_.countRejected(); return false;

  }

  // MQRANKSUM
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  mqrs_stats_.countUnfiltered();
  switch (checkInfoField(variant, MQRANKSUM_FIELD_, MQRANKSUM_LEVEL_, true)) {

    case FieldOutcome::MISSING: mqrs_stats_.countMissing();  break;
    case FieldOutcome::PASS:    mqrs_stats_.countAccepted(); break;
    case FieldOutcome::FAIL:    mqrs_stats_.countRejected(); return false;

  }

  // READPOSRANKSUM
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  rprs_stats_.countUnfiltered();
  switch (checkInfoField(variant, READPOSRANKSUM_FIELD_, READPOSRANKSUM_LEVEL_, true)) {

    case FieldOutcome::MISSING: rprs_stats_.countMissing();  break;
    case FieldOutcome::PASS:    rprs_stats_.countAccepted(); break;
    case FieldOutcome::FAIL:    rprs_stats_.countRejected(); return false;

  }

  // READDEPTH
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  if constexpr (READDEPTH_ACTIVE_) {

    depth_stats_.countUnfiltered();
    if (auto const& format_data = variant.evidence().formatData()) {

      if ((*format_data)->DPCount() >= MINIMUM_READDEPTH_) {

        depth_stats_.countAccepted();

      } else {

        depth_stats_.countRejected();
        return false;

      }

    } else {

      depth_stats_.countMissing();

    }

  } // if constexpr (READDEPTH_ACTIVE_)

  ++accepted_variants_;

  return true;

}
