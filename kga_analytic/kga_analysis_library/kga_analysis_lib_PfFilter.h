//
// Created by kellerberrin on 14/11/23.
//

#ifndef KGA_ANALYSIS_PFEMP_FILTER_H
#define KGA_ANALYSIS_PFEMP_FILTER_H


#include "kgl_pf7_sample_parser.h"
#include "kgl_pf7_fws_parser.h"
#include "kgl_pf3k_coi.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome::analysis {   //  organization::project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class FilterPf7 {

public:

  FilterPf7(std::shared_ptr<const Pf7SampleResource> Pf7_sample_ptr,
            std::shared_ptr<const Pf7FwsResource> Pf7_fws_ptr) : Pf7_sample_ptr_(std::move(Pf7_sample_ptr)),
                                                                 Pf7_fws_ptr_(std::move(Pf7_fws_ptr)) {}
  ~FilterPf7() = default;

  [[nodiscard]] std::shared_ptr<PopulationDB> qualityFilter( const std::shared_ptr<const PopulationDB>& unfiltered_population_ptr) const;

  // Get
  [[nodiscard]] bool codingFilterActive() const { return coding_filter_active_; } // Restrict to gene coding areas only.
  [[nodiscard]] bool filterQCActive() const { return filter_qc_active_; }
  [[nodiscard]] bool filterFWSActive() const { return filter_fws_active_; }
  [[nodiscard]] double FWSMonoclonalThreshold() const { return fws_monoclonal_threshold_; }
  [[nodiscard]] bool MleafFilterActive() const { return mleaf_filter_active_; }
  [[nodiscard]] bool AFFilterActive() const { return af_filter_active_; }
  [[nodiscard]] bool SNPFilterActive() const { return snp_filter_active_; }  // Only SNP variants
  [[nodiscard]] double VariantFrequencyCutoff() const { return variant_frequency_cutoff_; }
  // Set
  void codingFilterActive(bool coding_filter_active) { coding_filter_active_ = coding_filter_active; } // Restrict to gene coding areas only.
  void filterQCActive(bool filter_qc_active) { filter_qc_active_ = filter_qc_active; }
  void filterFWSActive(bool filter_fws_active) { filter_fws_active_ = filter_fws_active; }
  void FWSMonoclonalThreshold(double fws_monoclonal_threshold) { fws_monoclonal_threshold_ = fws_monoclonal_threshold; }
  void MleafFilterActive(bool mleaf_filter_active) { mleaf_filter_active_ = mleaf_filter_active; }
  void AFFilterActive(bool af_filter_active) { af_filter_active_ = af_filter_active; };
  void SNPFilterActive(bool snp_filter_active) { snp_filter_active_ = snp_filter_active; }  // Only SNP variants
  void VariantFrequencyCutoff(double variant_frequency_cutoff) { variant_frequency_cutoff_ = variant_frequency_cutoff; }

private:

  std::shared_ptr<const Pf7SampleResource> Pf7_sample_ptr_;
  std::shared_ptr<const Pf7FwsResource> Pf7_fws_ptr_;

  // Filter constants.
  bool coding_filter_active_{true}; // Restrict to gene coding areas only.
  bool filter_qc_active_{true};
  bool filter_fws_active_{true};
  double fws_monoclonal_threshold_{0.95};
  bool mleaf_filter_active_{false};
  bool af_filter_active_{false};
  bool snp_filter_active_{false};  // Only SNP variants
  double variant_frequency_cutoff_{0.0};


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FilterPf3k {

public:

  FilterPf3k(std::shared_ptr<const Pf3kCOIResource> Pf3k_COI_ptr) : Pf3k_COI_ptr_(std::move(Pf3k_COI_ptr)) {}
  ~FilterPf3k() = default;

  [[nodiscard]] std::shared_ptr<PopulationDB> filterCOI( const std::shared_ptr<const PopulationDB>& unfiltered_population_ptr) const;

private:

  std::shared_ptr<const Pf3kCOIResource> Pf3k_COI_ptr_;


};



} // Namespace


#endif //KGA_ANALYSIS_PFEMP_FILTER_H
