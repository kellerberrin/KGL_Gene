//
// Created by kellerberrin on 7/6/20.
//

#include "kgl_age_analysis.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"


namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///


double kgl::InfoAgeAnalysis::processField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(variant_ptr, field_name);

  if (field_opt) {

    std::vector<int64_t> field_vec = InfoEvidenceAnalysis::varianttoIntegers(field_opt.value());

    if (field_vec.size() != 1) {

      ExecEnv::log().warn("InfoAgeAnalysis::processVariant, Field: {} expected vector size 1, get vector size: {}",
                          field_name, field_vec.size());

    } else {

      return static_cast<double>(field_vec.front());

    }

  } else {

    ExecEnv::log().warn("InfoAgeAnalysis::processVariant, Field: {} not available", field_name);

  }

  return 0.0;

}



std::vector<double> kgl::InfoAgeAnalysis::processBin(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> bin_info_opt = InfoEvidenceAnalysis::getInfoData(variant_ptr, field_name);

  if (bin_info_opt) {

    std::vector<std::string> age_string = InfoEvidenceAnalysis::varianttoStrings(bin_info_opt.value());

    std::vector<double> age_vector = InfoEvidenceAnalysis::stringBinToFloat(age_string, AGE_BIN_SIZE_);

    return age_vector;

  } else {

    ExecEnv::log().warn("InfoAgeAnalysis::processVariant, data for age bin field: {}, not available", field_name);

  }

  return std::vector<double>(AGE_BIN_SIZE_, 0.0);

}


bool kgl::InfoAgeAnalysis::processVariant(const std::shared_ptr<const Variant>& variant_ptr) {

  ++variant_count_;

  het_under_30_ += processField(variant_ptr, HETERO_UNDER30_FIELD_);
  het_80_over_ += processField(variant_ptr, HETERO_80OVER_FIELD_);
  hom_under_30_ += processField(variant_ptr, HOMO_UNDER30_FIELD_);
  hom_80_over_ += processField(variant_ptr, HOMO_80OVER_FIELD_);

  size_t index = 0;
  for (auto& age : processBin(variant_ptr, HETERO_AGE_FIELD_)) {

    het_age_vector_[index] += age;

    ++index;

  }

  index = 0;
  for (auto& age : processBin(variant_ptr, HOMO_AGE_FIELD_)) {

    hom_age_vector_[index] += age;

    ++index;

  }

  return true;

}


void kgl::InfoAgeAnalysis::addAgeAnalysis(const InfoAgeAnalysis& age_analysis) {

  size_t index = 0;
  for (auto age : age_analysis.hom_age_vector_) {

    hom_age_vector_[index] += age;
    ++index;

  }

  hom_under_30_ += age_analysis.hom_under_30_;
  hom_80_over_ += age_analysis.hom_80_over_;

  index = 0;
  for (auto age : age_analysis.het_age_vector_) {

    het_age_vector_[index] += age;
    ++index;

  }

  het_under_30_ += age_analysis.het_under_30_;
  het_80_over_ += age_analysis.het_80_over_;

  variant_count_ += age_analysis.variant_count_;

}



double kgl::InfoAgeAnalysis::sumHomozygous() const {

  double sum = hom_under_30_;
  for (auto age : hom_age_vector_) {

    sum += age;

  }

  sum += hom_80_over_;

  return sum;

}

double kgl::InfoAgeAnalysis::sumHeterozygous() const {

  double sum = het_under_30_;
  for (auto age : het_age_vector_) {

    sum += age;

  }

  sum += het_80_over_;

  return sum;

}



// Utility function writes results to a stream.
std::ostream& operator<<(std::ostream& ostream, const kellerberrin::genome::InfoAgeAnalysis& age_analysis) {

  ostream << "Age Analysis, " << age_analysis.title() << ", Variant Count: " << age_analysis.variantCount() << '\n';

  ostream << kgl::InfoAgeAnalysis::header() << '\n';
  ostream << "hom, " << age_analysis.ageHomozygousUnder30() << ", ";
  for (auto const age : age_analysis.ageHomozygousVector()) {

    ostream << age << ", ";

  }
  ostream << age_analysis.ageHomozygous80Over() << '\n';

  double sum_hom = age_analysis.sumHomozygous() * 0.01;
  ostream << "hom%, " << (sum_hom > 0 ? (age_analysis.ageHomozygousUnder30() / sum_hom) : 0.0) << ", ";
  for (auto const age : age_analysis.ageHomozygousVector()) {

    ostream << (sum_hom > 0 ? (age /sum_hom) : 0.0) << ", ";

  }
  ostream << (sum_hom > 0 ? (age_analysis.ageHomozygous80Over() /sum_hom) : 0.0) << '\n';


  ostream << "het, " << age_analysis.ageHeterozygousUnder30() << ", ";
  for (auto const age : age_analysis.ageHeterozygousVector()) {

    ostream << age << ", ";

  }
  ostream << age_analysis.ageHeterozygous80Over() << '\n';

  double sum_het = age_analysis.sumHeterozygous() * 0.01;

  ostream << "het%, " << (sum_het > 0 ? (age_analysis.ageHeterozygousUnder30() / sum_het) : 0.0)  << ", ";
  for (auto const age : age_analysis.ageHeterozygousVector()) {

    ostream << (sum_het > 0 ? (age / sum_het) : 0.0) << ", ";

  }
  ostream << (sum_het > 0 ? (age_analysis.ageHeterozygous80Over() / sum_het) : 0.0) << '\n';

  if (age_analysis.ageHomozygousUnder30() > 0) {

    ostream << "het/hom, " << (age_analysis.ageHeterozygousUnder30() / age_analysis.ageHomozygousUnder30()) << ", ";

  } else {

    ostream << "het/hom, " << 0.0 << ", ";

  }
  size_t index = 0;
  for (auto const age : age_analysis.ageHomozygousVector()) {

    if (age > 0) {

      ostream << (age_analysis.ageHeterozygousVector()[index] / age) << ", ";

    } else {

      ostream << 0.0 << ", ";

    }

    ++index;

  }
  if (age_analysis.ageHomozygous80Over() > 0) {

    ostream << (age_analysis.ageHeterozygous80Over() / age_analysis.ageHomozygous80Over()) << '\n' << '\n';

  } else {

    ostream << 0.0 << '\n' << '\n';

  }


  return ostream;

}

