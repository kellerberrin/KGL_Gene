//
// Created by kellerberrin on 27/04/23.
//

#ifndef KGL_ANALYSIS_PFEMP_FILTER_H
#define KGL_ANALYSIS_PFEMP_FILTER_H


#include "kgl_variant.h"
#include "kgl_variant_db_population.h"
#include "kel_utility.h"
#include "kgl_variant_filter_type.h"


namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simple object to collect filter statistics for each filter field.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FilterFieldInfo {

public:

  explicit FilterFieldInfo(std::string filter_name) : filter_name_(std::move(filter_name)) {}
  ~FilterFieldInfo() = default;

  void initialize() { unfiltered_ = 0; accepted_ = 0; rejected_ = 0; missing_ = 0; }
  void printStats(double filter_level) const;


  std::atomic_uint64_t unfiltered_{0};
  std::atomic_uint64_t accepted_{0};
  std::atomic_uint64_t rejected_{0};
  std::atomic_uint64_t missing_{0};

private:

  const std::string filter_name_;

  [[nodiscard]] double acceptedPercent() const { return unfiltered_ > 0 ? 100.0 * (static_cast<double>(accepted_) / static_cast<double>(unfiltered_)) : 0.0; }
  [[nodiscard]] double rejectedPercent() const { return unfiltered_ > 0 ? 100.0 * (static_cast<double>(rejected_) / static_cast<double>(unfiltered_)) : 0.0; }
  [[nodiscard]] double missingPercent() const { return unfiltered_ > 0 ? 100.0 * (static_cast<double>(missing_) / static_cast<double>(unfiltered_)) : 0.0; }
  [[nodiscard]] double rejectionRate() const { return (accepted_ + rejected_) > 0 ? 100.0 * (static_cast<double>(rejected_) / static_cast<double>((accepted_ + rejected_))) : 0.0; }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Frequency filter.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Which field to use.
enum class FreqInfoField { AF, MLEAF };

class P7FrequencyFilter : public FilterVariants {

public:

  P7FrequencyFilter(FreqInfoField info_field, double freq_cutoff) : info_field_(info_field), freq_cutoff_(freq_cutoff) { filterName("P7FrequencyFilter"); }
  ~P7FrequencyFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant &variant) const override;

  [[nodiscard]] std::shared_ptr <BaseFilter> clone() const override { return std::make_shared<P7FrequencyFilter>(*this); }

  constexpr static const char* INFO_FREQ_FILTER_{"InfoFreqFilter"};

  static void initializeStats() { freq_filter_stats_.initialize(); }
  static void printStats(double freq_cutoff) { freq_filter_stats_.printStats(freq_cutoff); }

private:

  constexpr static const char* AF_FIELD_{"AF"};
  constexpr static const char* MLEAF_FIELD_{"MLEAF"};

  FreqInfoField info_field_;
  double freq_cutoff_;

  inline static FilterFieldInfo freq_filter_stats_{INFO_FREQ_FILTER_};

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Bespoke variant quality filter for the P7 Pf database.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class P7VariantFilter : public FilterVariants {

public:

  P7VariantFilter() { filterName("P7VariantFilter"); }
  ~P7VariantFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant &variant) const override;

  [[nodiscard]] std::shared_ptr <BaseFilter> clone() const override { return std::make_shared<P7VariantFilter>(*this); }

  // Before filtering.
  static void initializeStats();
  // After filtering.
  static void printStats();

private:

  // Filter constants.
  constexpr static const char* VQSLOD_FIELD_{"VQSLOD"};
  // constexpr static const double VQSLOD_LEVEL_{1.775};
  // constexpr static const double VQSLOD_LEVEL_{1.2168};
  constexpr static const double VQSLOD_LEVEL_{0.0};
  constexpr static const char* QD_FIELD_{"QD"};
  constexpr static const double QD_LEVEL_{2.0};
  constexpr static const char* MQ_FIELD_{"MQ"};
  constexpr static const double MQ_LEVEL_{30.0};
  constexpr static const char* Pf7_Pf3D7_MITO{"Pf3D7_MIT_v3"};
  constexpr static const double Pf7_Pf3D7_MITO_MQ_LEVEL_{5.0};
  constexpr static const char* SOR_FIELD_{"SOR"};
  constexpr static const double SOR_LEVEL_{3.0};
  constexpr static const char* MQRANKSUM_FIELD_{"MQRankSum"};
  constexpr static const double MQRANKSUM_LEVEL_{-12.5};
  constexpr static const char* READPOSRANKSUM_FIELD_{"ReadPosRankSum"};
  constexpr static const double READPOSRANKSUM_LEVEL_{-8.0};
  constexpr static const char* READDEPTH_FIELD_{"ReadDepth"};
  constexpr static const size_t MINIMUM_READDEPTH_{5};
  constexpr static const bool READDEPTH_ACTIVE_{false};


  // Global filter statistics.
  static inline std::atomic_uint64_t unfiltered_variants_{0};
  static inline std::atomic_uint64_t accepted_variants_{0};

  // The statistics fields are static with atomic counters.
  // This is necessary because filters are cloned and run in separate threads,
  inline static FilterFieldInfo vqslod_stats_{VQSLOD_FIELD_};
  inline static FilterFieldInfo qd_stats_{QD_FIELD_};
  inline static FilterFieldInfo mq_stats_{MQ_FIELD_};
  inline static FilterFieldInfo sor_stats_{SOR_FIELD_};
  inline static FilterFieldInfo mqrs_stats_{MQRANKSUM_FIELD_};
  inline static FilterFieldInfo rprs_stats_{READPOSRANKSUM_FIELD_};
  inline static FilterFieldInfo depth_stats_{READDEPTH_FIELD_};

};


} // Namespace.


#endif //KGL_ANALYSIS_PFEMP_FILTER_H
