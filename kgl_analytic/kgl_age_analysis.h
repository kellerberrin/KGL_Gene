//
// Created by kellerberrin on 7/6/20.
//

#ifndef KGL_AGE_ANALYSIS_H
#define KGL_AGE_ANALYSIS_H

#include "kgl_variant.h"


namespace kellerberrin::genome {   //  organization::project level namespace

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

  [[nodiscard]] double ageHomozygousUnder30() const { return hom_under_30_; }
  [[nodiscard]] double ageHomozygous80Over() const { return hom_80_over_; }
  [[nodiscard]] const std::vector<double>& ageHomozygousVector() const { return hom_age_vector_; }
  [[nodiscard]] const std::vector<double>& ageHeterozygousVector() const { return het_age_vector_; }
  [[nodiscard]] double ageHeterozygousUnder30() const { return het_under_30_; }
  [[nodiscard]] double ageHeterozygous80Over() const { return het_80_over_; }
  [[nodiscard]] double sumHomozygous() const;
  [[nodiscard]] double sumHeterozygous() const;

  [[nodiscard]] std::string title() const { return analysis_title_; }

  [[nodiscard]] std::string static header() { return "hom/het, <30, >=30, >=35, >=40, >=45, >=50, >=55, >=60, >=65, >=70, >=75, >=80"; }

  [[nodiscard]] size_t variantCount() const { return variant_count_; }

private:

  std::string analysis_title_;

  std::vector<double> hom_age_vector_;
  double hom_under_30_{0};
  double hom_80_over_{0};

  std::vector<double> het_age_vector_;
  double het_under_30_{0};
  double het_80_over_{0};

  size_t variant_count_{0};

  // Heterozygous bin field
  constexpr static const char* HETERO_AGE_FIELD_{"age_hist_het_bin_freq"};
  constexpr static const char* HETERO_UNDER30_FIELD_{"age_hist_het_n_smaller"};
  constexpr static const char* HETERO_80OVER_FIELD_{"age_hist_het_n_larger"};

  // Homozygous bin field
  constexpr static const char* HOMO_AGE_FIELD_{"age_hist_hom_bin_freq"};
  constexpr static const char* HOMO_UNDER30_FIELD_{"age_hist_hom_n_smaller"};
  constexpr static const char* HOMO_80OVER_FIELD_{"age_hist_hom_n_larger"};
  // Bin field size.
  constexpr static const size_t AGE_BIN_SIZE_{10};

  double processField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name);
  std::vector<double> processBin(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name);

};


} // namespace


std::ostream& operator<<(std::ostream& ostream, const kellerberrin::genome::InfoAgeAnalysis& age_analysis);



#endif //KGL_KGL_AGE_ANALYSIS_H
