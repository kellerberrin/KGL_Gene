//
// Created by kellerberrin on 7/6/20.
//

#ifndef KGL_AGE_ANALYSIS_H
#define KGL_AGE_ANALYSIS_H

#include "kgl_variant.h"


namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A matrix of heterozygous vs homozygous variant counts

class HetHomAgeMatrix {

public:

  HetHomAgeMatrix() {

    std::vector<double> row_vector(AGE_MATRIX_SIZE_, 0.0);

    for (size_t index = 0; index < AGE_MATRIX_SIZE_; ++index) {

      het_hom_age_matrix_.push_back(row_vector);

    }

  }
  ~HetHomAgeMatrix() = default;

  void updateMatrix(size_t het_age_index, size_t hom_age_index, double value) { het_hom_age_matrix_[het_age_index][hom_age_index] += value; }
  [[nodiscard]] const std::vector<std::vector<double>>& getMatrix() const { return het_hom_age_matrix_; }

private:

  constexpr static const size_t AGE_MATRIX_SIZE_{12};

  std::vector<std::vector<double>> het_hom_age_matrix_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to harvest age bin statistics (if they exist)

class InfoAgeAnalysis {

public:

  InfoAgeAnalysis(std::string analysis_title) : analysis_title_(std::move(analysis_title)),
                                                hom_age_vector_(AGE_BIN_SIZE_, 0),
                                                het_age_vector_(AGE_BIN_SIZE_, 0) {}
  InfoAgeAnalysis(const InfoAgeAnalysis&) = default;
  ~InfoAgeAnalysis() = default;

  void addAgeAnalysis(const InfoAgeAnalysis& age_analysis);

  bool processVariant(const std::shared_ptr<const Variant>& variant_ptr);

  [[nodiscard]] const std::vector<double>& ageHomozygousVector() const { return hom_age_vector_; }
  [[nodiscard]] const std::vector<double>& ageHeterozygousVector() const { return het_age_vector_; }
  [[nodiscard]] double sumHomozygous() const;
  [[nodiscard]] double sumHeterozygous() const;
  [[nodiscard]] double heteroHomoRatioAll() const { return sumHomozygous() > 0 ? sumHeterozygous() / sumHomozygous() : 0.0; }
  [[nodiscard]] double averageHomozygousAge() const { return sumHomozygous() > 0 ? (ageWeightedSumHomozygous() / sumHomozygous()) : 0.0; }
  [[nodiscard]] double averageHeterozygousAge() const { return sumHeterozygous() > 0 ? (ageWeightedSumHeterozygous() / sumHeterozygous()) : 0.0; }
  [[nodiscard]] double averageCombinedAge() const {
    return (sumHeterozygous() + sumHomozygous()) > 0 ? (ageWeightedSumHomozygous() + ageWeightedSumHeterozygous()) / (sumHeterozygous() + sumHomozygous()) : 0.0; }

  [[nodiscard]] std::string title() const { return analysis_title_; }

  [[nodiscard]] std::string static header() { return "hom/het, <30, >=30, >=35, >=40, >=45, >=50, >=55, >=60, >=65, >=70, >=75, >=80"; }

  [[nodiscard]] size_t variantCount() const { return variant_count_; }

  [[nodiscard]] const HetHomAgeMatrix& hetHomMatrix() const { return het_hom_matrix_; }

private:

  std::string analysis_title_;

  std::vector<double> hom_age_vector_;
  std::vector<double> het_age_vector_;

  double all_allele_{0.0};
  double all_alternate_allele_{0.0};

  size_t variant_count_{0};

  HetHomAgeMatrix het_hom_matrix_;

  // Heterozygous bin field
  constexpr static const char* HETERO_AGE_FIELD_{"age_hist_het_bin_freq"};
  constexpr static const char* HETERO_UNDER30_FIELD_{"age_hist_het_n_smaller"};
  constexpr static const char* HETERO_80OVER_FIELD_{"age_hist_het_n_larger"};

  // Homozygous bin field
  constexpr static const char* HOMO_AGE_FIELD_{"age_hist_hom_bin_freq"};
  constexpr static const char* HOMO_UNDER30_FIELD_{"age_hist_hom_n_smaller"};
  constexpr static const char* HOMO_80OVER_FIELD_{"age_hist_hom_n_larger"};

  constexpr static const char* TOTAL_ALLELE_COUNT_{"AN"};  // Total number of samples including reference alleles.
  constexpr static const char* ALTERNATE_ALLELE_COUNT_{"AC"}; // Total number of samples with alternate alleles.
  // Bin field size.
  constexpr static const size_t FIELD_AGE_BIN_SIZE_{10};
  constexpr static const size_t AGE_BIN_SIZE_{12};

  // Age weighting vector is the mid-point between the age bins.
  // Which appears to be a reasonable assumption.
  const std::vector<double> age_weight_vector_{ 15.0, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 85.0};

  double processField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name);
  std::vector<double> processBin(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name);

  [[nodiscard]] double ageWeightedSumHomozygous() const;
  [[nodiscard]] double ageWeightedSumHeterozygous() const;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Average Age sorted multimap of variants.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using AgeMultiMap = std::multimap<double, std::shared_ptr<const InfoAgeAnalysis>>;
class AgeSortedMap {

public:

  AgeSortedMap(AgeMultiMap&& homozygous_map, AgeMultiMap&& heterozygous_map, AgeMultiMap&& combined_map)
  : heterozygous_map_(heterozygous_map), homozygous_map_(homozygous_map), combined_map_(combined_map) {}
  AgeSortedMap() = default;
  ~AgeSortedMap() = default;

  [[nodiscard]] bool processVariant(const std::shared_ptr<const Variant>& variant_ptr);

  [[nodiscard]] std::pair<double, double> heterozygousAgeStatistics() const { return ageStatistics(heterozygous_map_); }
  [[nodiscard]] std::pair<double, double> homozygousAgeStatistics() const  { return ageStatistics(homozygous_map_); }
  [[nodiscard]] std::pair<double, double> combinedAgeStatistics() const  { return ageStatistics(combined_map_); }

  // Calling with (1.0, 0.9) or (100%, 90%) retrieves the 10% of variants with the highest average age.
  // Calling with (0.1, 0.0) or (10%, 0%) retrieves the 10% of variants with the lowest average age.
  [[nodiscard]] AgeMultiMap ageFilterHeterozygous(double top_percentile, double bottom_percentile) const
      { return ageFilter(top_percentile, bottom_percentile, heterozygous_map_); }
  [[nodiscard]] AgeMultiMap ageFilterHomozygous(double top_percentile, double bottom_percentile) const
      { return ageFilter(top_percentile, bottom_percentile, homozygous_map_); }
  [[nodiscard]] AgeMultiMap ageFilterCombined(double top_percentile, double bottom_percentile) const
      { return ageFilter(top_percentile, bottom_percentile, combined_map_); }

  [[nodiscard]] InfoAgeAnalysis aggregateHeteroAnalysis() const { return aggregateAnalysis(heterozygous_map_, "Heterozygous Aggregation"); }
  [[nodiscard]] InfoAgeAnalysis aggregateHomoAnalysis() const  { return aggregateAnalysis(homozygous_map_, "Homozygous Aggregation"); }
  [[nodiscard]] InfoAgeAnalysis aggregateCombinedAnalysis() const { return aggregateAnalysis(combined_map_, "Combined Aggregation"); }

private:

  AgeMultiMap heterozygous_map_;
  AgeMultiMap homozygous_map_;
  AgeMultiMap combined_map_;

  [[nodiscard]] AgeMultiMap ageFilter(double top_percentile, double bottom_percentile, const AgeMultiMap& age_map) const;
  [[nodiscard]] std::pair<double, double> ageStatistics(const AgeMultiMap& age_map) const;
  [[nodiscard]] InfoAgeAnalysis aggregateAnalysis(const AgeMultiMap& age_map, const std::string& title) const;

};




} // namespace


std::ostream& operator<<(std::ostream& ostream, const kellerberrin::genome::InfoAgeAnalysis& age_analysis);

std::ostream& operator<<(std::ostream& ostream, const kellerberrin::genome::AgeSortedMap& sorted_age_analysis);



#endif //KGL_KGL_AGE_ANALYSIS_H
