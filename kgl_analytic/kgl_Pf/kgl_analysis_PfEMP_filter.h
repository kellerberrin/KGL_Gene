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

  FilterFieldInfo(std::string filter_name, double filter_level) : filter_name_(std::move(filter_name)), filter_level_(filter_level) {}
  ~FilterFieldInfo() = default;

  void initialize() { unfiltered_ = 0; accepted_ = 0; rejected_ = 0; missing_ = 0; }
  void printStats() const;


  std::atomic_uint64_t unfiltered_{0};
  std::atomic_uint64_t accepted_{0};
  std::atomic_uint64_t rejected_{0};
  std::atomic_uint64_t missing_{0};

private:

  const std::string filter_name_;
  double filter_level_{0.0};

  [[nodiscard]] double acceptedPercent() const { return unfiltered_ > 0 ? 100.0 * (static_cast<double>(accepted_) / static_cast<double>(unfiltered_)) : 0.0; }
  [[nodiscard]] double rejectedPercent() const { return unfiltered_ > 0 ? 100.0 * (static_cast<double>(rejected_) / static_cast<double>(unfiltered_)) : 0.0; }
  [[nodiscard]] double missingPercent() const { return unfiltered_ > 0 ? 100.0 * (static_cast<double>(missing_) / static_cast<double>(unfiltered_)) : 0.0; }
  [[nodiscard]] double rejectionRate() const { return (accepted_ + rejected_) > 0 ? 100.0 * (static_cast<double>(rejected_) / static_cast<double>((accepted_ + rejected_))) : 0.0; }

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
  constexpr static const double VQSLOD_LEVEL_{0.5};
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
  inline static FilterFieldInfo vqslod_stats_{VQSLOD_FIELD_, VQSLOD_LEVEL_};
  inline static FilterFieldInfo qd_stats_{QD_FIELD_, QD_LEVEL_};
  inline static FilterFieldInfo mq_stats_{MQ_FIELD_, MQ_LEVEL_};
  inline static FilterFieldInfo sor_stats_{SOR_FIELD_, SOR_LEVEL_};
  inline static FilterFieldInfo mqrs_stats_{MQRANKSUM_FIELD_, MQRANKSUM_LEVEL_};
  inline static FilterFieldInfo rprs_stats_{READPOSRANKSUM_FIELD_, READPOSRANKSUM_LEVEL_};
  inline static FilterFieldInfo depth_stats_{READDEPTH_FIELD_, MINIMUM_READDEPTH_};

};


} // Namespace.


#endif //KGL_ANALYSIS_PFEMP_FILTER_H
