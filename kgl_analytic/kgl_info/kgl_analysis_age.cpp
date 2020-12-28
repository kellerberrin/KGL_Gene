//
// Created by kellerberrin on 7/6/20.
//

#include "kgl_analysis_age.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"

#include <numeric>

namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///


double kgl::InfoAgeAnalysis::processField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, field_name);

  if (field_opt) {

    std::vector<int64_t> field_vec = InfoEvidenceAnalysis::varianttoIntegers(field_opt.value());

    if (field_vec.size() != 1) {

      ExecEnv::log().warn("InfoAgeAnalysis::processVariant, Field: {} expected vector size 1, get vector size: {}",
                          field_name, field_vec.size());
      return 0.0;

    } else {

      return static_cast<double>(field_vec.front());

    }

  }

  return 0.0;

}



std::vector<double> kgl::InfoAgeAnalysis::processBin(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> bin_info_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, field_name);

  if (bin_info_opt) {

    std::vector<std::string> age_string = InfoEvidenceAnalysis::varianttoStrings(bin_info_opt.value());

    std::vector<double> age_vector = InfoEvidenceAnalysis::stringBinToFloat(age_string, FIELD_AGE_BIN_SIZE_);

    return age_vector;

  }

//  ExecEnv::log().warn("InfoAgeAnalysis::processVariant, data for age bin field: {}, not available", field_name);

  return std::vector<double>(FIELD_AGE_BIN_SIZE_, 0.0);

}


bool kgl::InfoAgeAnalysis::processVariant(const std::shared_ptr<const Variant>& variant_ptr) {

  ++variant_count_;

  double het_under_30 = processField(variant_ptr, HETERO_UNDER30_FIELD_);
  double het_80_over = processField(variant_ptr, HETERO_80OVER_FIELD_);
  auto var_het_vector = processBin(variant_ptr, HETERO_AGE_FIELD_);

  std::vector<double> het_age_vector;
  het_age_vector.push_back(het_under_30);
  for (auto age : var_het_vector) {

    het_age_vector.push_back(age);

  }
  het_age_vector.push_back(het_80_over );

  if (het_age_vector.size() != AGE_BIN_SIZE_) {

    ExecEnv::log().error("InfoAgeAnalysis::processVariant, het_age_vector size: {}, expected size: {}", het_age_vector.size(), AGE_BIN_SIZE_);
    return false;

  }

  double hom_under_30 = processField(variant_ptr, HOMO_UNDER30_FIELD_);
  double hom_80_over = processField(variant_ptr, HOMO_80OVER_FIELD_);
  auto var_hom_vector = processBin(variant_ptr, HOMO_AGE_FIELD_);

  std::vector<double> hom_age_vector;
  hom_age_vector.push_back(hom_under_30);
  for (auto age : var_hom_vector) {

    hom_age_vector.push_back(age);

  }
  hom_age_vector.push_back(hom_80_over);

  if (hom_age_vector.size() != AGE_BIN_SIZE_) {

    ExecEnv::log().error("InfoAgeAnalysis::processVariant, hom_age_vector size: {}, expected size: {}", hom_age_vector.size(), AGE_BIN_SIZE_);
    return false;

  }

  std::transform (het_age_vector_.begin(),
             het_age_vector_.end(),
                   het_age_vector.begin(),
                   het_age_vector_.begin(),
                   std::plus<double>());

  std::transform (hom_age_vector_.begin(),
                   hom_age_vector_.end(),
                   hom_age_vector.begin(),
                   hom_age_vector_.begin(),
                   std::plus<double>());

  all_allele_ += processField(variant_ptr, TOTAL_ALLELE_COUNT_);
  all_alternate_allele_ += processField(variant_ptr, ALTERNATE_ALLELE_COUNT_);

// normalize the vector.
  auto het_sum = std::accumulate(het_age_vector.begin(), het_age_vector.end(), decltype(het_age_vector)::value_type(0.0));
  for (size_t index = 0; index < AGE_BIN_SIZE_; ++index) {

    if (het_sum > 0) {

      het_age_vector[index] = het_age_vector[index] / het_sum;

    } else {

      het_age_vector[index] = 0.0;

    }

  }

  double check_sum = 0.0;
  for (auto const& age: het_age_vector) {

    check_sum += age;

  }

  if (check_sum != 0.0 and (check_sum < 0.999999 or check_sum > 1.000001)) {

    ExecEnv::log().error("InfoAgeAnalysis::processVariant, Normalized sum(het_age_vector) out of bounds: {}, expected 1.0", check_sum);

  }

  auto hom_sum = std::accumulate(hom_age_vector.begin(), hom_age_vector.end(), decltype(hom_age_vector)::value_type(0.0));
  for (size_t index = 0; index < AGE_BIN_SIZE_; ++index) {

    if (hom_sum > 0) {

      hom_age_vector[index] = hom_age_vector[index] / hom_sum;

    } else {

      hom_age_vector[index] = 0.0;

    }

  }

  check_sum = 0.0;
  for (auto const& age: hom_age_vector) {

    check_sum += age;

  }

  if (check_sum != 0.0 and (check_sum < 0.999999 or check_sum > 1.000001)) {

    ExecEnv::log().error("InfoAgeAnalysis::processVariant, Normalized sum(hom_age_vector) out of bounds: {}, expected 1.0", check_sum);

  }

  // Update the het/hom matrix
  check_sum = 0.0;
  for (size_t het_index = 0; het_index < AGE_BIN_SIZE_; ++het_index) {

    for (size_t hom_index = 0; hom_index < AGE_BIN_SIZE_; ++hom_index) {

      double cross_product = het_age_vector[het_index] * hom_age_vector[hom_index];
      het_hom_matrix_.updateMatrix(het_index, hom_index, cross_product);
      check_sum += cross_product;

    }

  }

  if (check_sum != 0.0 and (check_sum < 0.999999 or check_sum > 1.000001)) {

    ExecEnv::log().error("InfoAgeAnalysis::processVariant, Normalized sum(het_hom_matrix) out of bounds: {}, expected 1.0", check_sum);

  }

  return true;

}


void kgl::InfoAgeAnalysis::addAgeAnalysis(const InfoAgeAnalysis& age_analysis) {

  size_t index = 0;
  for (auto age : age_analysis.hom_age_vector_) {

    hom_age_vector_[index] += age;
    ++index;

  }

  index = 0;
  for (auto age : age_analysis.het_age_vector_) {

    het_age_vector_[index] += age;
    ++index;

  }

  variant_count_ += age_analysis.variant_count_;

}


double kgl::InfoAgeAnalysis::sumHomozygous() const {

  double sum = 0;
  for (auto age : hom_age_vector_) {

    sum += age;

  }

  return sum;

}

double kgl::InfoAgeAnalysis::sumHeterozygous() const {

  double sum = 0;
  for (auto age : het_age_vector_) {

    sum += age;

  }

  return sum;

}


double kgl::InfoAgeAnalysis::ageWeightedSumHomozygous() const {

  size_t index = 0;
  double sum = 0.0;
  for (auto age : hom_age_vector_) {

    sum += age * age_weight_vector_[index];
    ++index;

  }

  return sum;

}

double kgl::InfoAgeAnalysis::ageWeightedSumHeterozygous() const {

  double sum = 0;
  size_t index = 0;
  for (auto age : het_age_vector_) {

    sum += age * age_weight_vector_[index];
    ++index;

  }

  return sum;

}



// Utility function writes results to a stream.
std::ostream& operator<<(std::ostream& ostream, const kellerberrin::genome::InfoAgeAnalysis& age_analysis) {

  ostream << "Age Analysis, " << age_analysis.title() << ", Variant Count: " << age_analysis.variantCount()
          << ", het/hom ratio all: " << age_analysis.heteroHomoRatioAll() << '\n';
  ostream << "hom Av Age: " <<  age_analysis.averageHomozygousAge() << ", het Av Age: " << age_analysis.averageHeterozygousAge() << '\n';

  ostream << kgl::InfoAgeAnalysis::header() << '\n';
  ostream << "hom, ";
  for (auto const age : age_analysis.ageHomozygousVector()) {

    ostream << age << ", ";

  }

  ostream << '\n';

  double sum_hom = age_analysis.sumHomozygous() * 0.01;
  ostream << "hom%,";
  for (auto const age : age_analysis.ageHomozygousVector()) {

    ostream << (sum_hom > 0 ? (age /sum_hom) : 0.0) << ", ";

  }

  ostream << '\n';

  ostream << "het, ";
  for (auto const age : age_analysis.ageHeterozygousVector()) {

    ostream << age << ", ";

  }

  ostream << '\n';


  ostream << "het%, ";
  double sum_het = age_analysis.sumHeterozygous() * 0.01;
  for (auto const age : age_analysis.ageHeterozygousVector()) {

    ostream << (sum_het > 0 ? (age / sum_het) : 0.0) << ", ";

  }

  ostream << '\n';

  ostream << "het/hom, ";
  size_t index = 0;
  for (auto const age : age_analysis.ageHomozygousVector()) {

    if (age > 0) {

      ostream << (age_analysis.ageHeterozygousVector()[index] / age) << ", ";

    } else {

      ostream << 0.0 << ", ";

    }

    ++index;

  }

  ostream << '\n' << '\n';
  ostream << kgl::InfoAgeAnalysis::header() << '\n';

  double matrix_sum = age_analysis.hetHomMatrix().sum();
  // Update the het/hom matrix
  for (size_t het_index = 0; het_index < age_analysis.hetHomMatrix().getMatrix().size(); ++het_index) {

    for (size_t hom_index = 0; hom_index < age_analysis.hetHomMatrix().getMatrix().size(); ++hom_index) {

      if (matrix_sum > 0) {

        ostream << (age_analysis.hetHomMatrix().getMatrix()[het_index][hom_index] / matrix_sum) * 100.0 << ',';

      } else {

        ostream << age_analysis.hetHomMatrix().getMatrix()[het_index][hom_index] << ',';

      }

    }

    ostream << '\n';

  }

  ostream << '\n' << std::endl;

  return ostream;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Average Age sorted multimap of variants.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::AgeSortedMap::processVariant(const std::shared_ptr<const Variant>& variant_ptr) {

  std::shared_ptr<InfoAgeAnalysis> age_analysis_ptr(std::make_shared<InfoAgeAnalysis>("Age Sorted Analysis"));

  age_analysis_ptr->processVariant(variant_ptr);

  if (age_analysis_ptr->averageHomozygousAge() > 0.0) {

    homozygous_map_.emplace(age_analysis_ptr->averageHomozygousAge(), age_analysis_ptr);

  }

  if (age_analysis_ptr->averageHeterozygousAge() > 0.0) {

    heterozygous_map_.emplace(age_analysis_ptr->averageHeterozygousAge(), age_analysis_ptr);

  }

  if (age_analysis_ptr->averageCombinedAge() > 0.0) {

    combined_map_.emplace(age_analysis_ptr->averageCombinedAge(), age_analysis_ptr);

  }

  return true;

}


kgl::AgeMultiMap kgl::AgeSortedMap::ageFilter(double top_percentile, double bottom_percentile, const AgeMultiMap& age_map) const {

  auto size = static_cast<double>(age_map.size());

  auto top_count = static_cast<size_t>(size * top_percentile);
  auto bottom_count = static_cast<size_t>(size * bottom_percentile);

  size_t count = 0;
  AgeMultiMap filtered_map;
  for (auto const& [average_age, age_analysis_ptr] :  age_map) {

    if (count >= bottom_count and count <= top_count) {

      filtered_map.emplace(average_age, age_analysis_ptr);

    }

    ++count;

  }

  return filtered_map;

}


std::pair<double, double> kgl::AgeSortedMap::ageStatistics(const AgeMultiMap& age_map) const {

  std::vector<double> average_age_vector;
  average_age_vector.reserve(age_map.size());

  for (auto const& age_analysis : age_map) {

    average_age_vector.push_back(age_analysis.first);

  }

  return Utility::stddev(average_age_vector);

}


kgl::InfoAgeAnalysis kgl::AgeSortedMap::aggregateAnalysis(const AgeMultiMap& age_map,const std::string& title) const {

  InfoAgeAnalysis aggregate_analysis(title);

  // Only need to aggregate one map as both maps contain
  // the age analysis for all variants.
  for (auto const& age_analysis : age_map) {

    aggregate_analysis.addAgeAnalysis(*age_analysis.second);

  }

  return aggregate_analysis;

}



// Formatted output to file.
std::ostream& operator<<(std::ostream& ostream, const kellerberrin::genome::AgeSortedMap& age_sorted_map) {


  std::pair<double, double> hom_stats = age_sorted_map.homozygousAgeStatistics();
  std::pair<double, double> het_stats = age_sorted_map.heterozygousAgeStatistics();
  std::pair<double, double> comb_stats = age_sorted_map.combinedAgeStatistics();
  ostream << "All ages Homozygous mean: " << hom_stats.first << ", stdev: " << hom_stats.second << '\n';
  ostream << "All ages Heterozygous mean: " << het_stats.first << ", stdev: " << het_stats.second << '\n';
  ostream << "All ages Combined mean: " << comb_stats.first << ", stdev: " << comb_stats.second << '\n';
  ostream << "All ages Homozygous aggregate results" << '\n';
  ostream << age_sorted_map.aggregateHomoAnalysis();
  ostream << "All ages Heterozygous aggregate results" << '\n';
  ostream << age_sorted_map.aggregateHeteroAnalysis();
  ostream << "All ages Combined aggregate results" << '\n';
  ostream << age_sorted_map.aggregateCombinedAnalysis();

  // Filter into age deciles and report results.
  for (double lower_filter = 0.0; lower_filter < 0.99; lower_filter += 0.1) {

    double upper_filter = lower_filter + 0.1;
    kgl::AgeSortedMap filtered_age_sorted_map( age_sorted_map.ageFilterHomozygous(upper_filter, lower_filter),
                                               age_sorted_map.ageFilterHeterozygous(upper_filter, lower_filter),
                                               age_sorted_map.ageFilterCombined(upper_filter, lower_filter));
    hom_stats = filtered_age_sorted_map.homozygousAgeStatistics();
    ostream << "Filter from :" << lower_filter * 100 << "%, to " << upper_filter * 100 << "%, Homozygous mean: "
            << hom_stats.first << ", stdev: " << hom_stats.second << '\n';
    het_stats = filtered_age_sorted_map.heterozygousAgeStatistics();
    ostream << "Filter from :" << lower_filter * 100 << "%, to " << upper_filter * 100 << "%, Heterozygous mean: "
            << het_stats.first << ", stdev: " << het_stats.second << '\n';
    comb_stats = filtered_age_sorted_map.combinedAgeStatistics();
    ostream << "Filter from :" << lower_filter * 100 << "%, to " << upper_filter * 100 << "%, Combined mean: "
            << comb_stats.first << ", stdev: " << comb_stats.second << '\n';
    ostream << "Filter from :" << lower_filter * 100 << "%, to " << upper_filter * 100 << "%, Homozygous Aggregate Analysis" << '\n';
    ostream << filtered_age_sorted_map.aggregateHomoAnalysis();
    ostream << "Filter from :" << lower_filter * 100 << "%, to " << upper_filter * 100 << "%, Heterozygous Aggregate Analysis" << '\n';
    ostream << filtered_age_sorted_map.aggregateHeteroAnalysis();
    ostream << "Filter from :" << lower_filter * 100 << "%, to " << upper_filter * 100 << "%, Combined Aggregate Analysis" << '\n';
    ostream << filtered_age_sorted_map.aggregateCombinedAnalysis();

  }

  // flush the results.
  ostream << std::endl;

  return ostream;

}
