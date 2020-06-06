//
// Created by kellerberrin on 4/6/20.
//

#ifndef KGL_ANALYSIS_INFO_FILTER_H
#define KGL_ANALYSIS_INFO_FILTER_H


#include "kgl_runtime.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_unphased_population.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"

#include "kgl_analysis_null.h"
#include "kgl_filter.h"


namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to harvest age bin statistics (if they exist)

class InfoAgeAnalysis {

public:

  InfoAgeAnalysis() : hom_age_vector_(AGE_BIN_SIZE_, 0), het_age_vector_(AGE_BIN_SIZE_, 0) {}
  ~InfoAgeAnalysis() = default;

  void addAgeAnalysis(const InfoAgeAnalysis& age_analysis);

  void processVariant(const std::shared_ptr<const Variant>& variant_ptr);

  [[nodiscard]] double ageHomozygousUnder30() const { return hom_under_30_; }
  [[nodiscard]] double ageHomozygous80Over() const { return hom_80_over_; }
  [[nodiscard]] const std::vector<double>& ageHomozygousVector() const { return hom_age_vector_; }
  [[nodiscard]] const std::vector<double>& ageHeterozygousVector() const { return het_age_vector_; }
  [[nodiscard]] double ageHeterozygousUnder30() const { return het_under_30_; }
  [[nodiscard]] double ageHeterozygous80Over() const { return het_80_over_; }
  [[nodiscard]] double sumHomozygous() const;
  [[nodiscard]] double sumHeterozygous() const;

  [[nodiscard]] std::string static header() { return "hom/het, <30, >=30, >=35, >=40, >=45, >=50, >= 55, >=60, >=65, >=70, >=75, >=80"; }

  [[nodiscard]] size_t variantCount() const { return variant_count_; }

private:

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



// Test the Info data filters.
class InfoFilterAnalysis : public NullAnalysis {

public:

  InfoFilterAnalysis() = default;
  ~InfoFilterAnalysis() = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "INFO_FILTER"; }
  [[nodiscard]] std::unique_ptr<NullAnalysis> factory() const override { return std::make_unique<InfoFilterAnalysis>(); }

  // Setup the analytics to process VCF data.
  // This function must be redefined.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const RuntimeParameterMap& named_parameters,
                                         std::shared_ptr<const GenomeCollection> reference_genomes) override;

  // Perform the genetic analysis per VCF file (need not be redefined)
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const UnphasedPopulation> vcf_population) override;

  // Perform the genetic analysis per iteration (need not be redefined)
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results (need not be redefined)
  [[nodiscard]] bool finalizeAnalysis() override;

  [[nodiscard]] bool getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters);

private:

  InfoAgeAnalysis age_analysis_all_;
  InfoAgeAnalysis age_analysis_0_1_percent_;
  InfoAgeAnalysis age_analysis_1_percent_;
  InfoAgeAnalysis age_analysis_2_percent_;
  InfoAgeAnalysis age_analysis_5_percent_;
  InfoAgeAnalysis age_analysis_10_percent_;

  std::string output_file_name_;

  constexpr static const char* OUTPUT_FILE_ = "OUTPUTFILE";

};



} // namespace


std::ostream& operator<<(std::ostream& ostream, const kellerberrin::genome::InfoAgeAnalysis& age_analysis);


#endif //KGL_ANALYSIS_INFO_FILTER_H
