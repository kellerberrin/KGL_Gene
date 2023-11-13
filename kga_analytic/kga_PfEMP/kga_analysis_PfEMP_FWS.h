//
// Created by kellerberrin on 22/05/23.
//

#ifndef KGL_ANALYSIS_PFEMP_FWS_H
#define KGL_ANALYSIS_PFEMP_FWS_H


#include "kga_analysis_PfEMP_heterozygous.h"
#include "kgl_variant_db_variant.h"


namespace kellerberrin::genome::analysis {   //  organization::project level namespace

enum class AlleleFrequencyBins: size_t { PERCENT_0_5 = 0,
                                PERCENT_5_10 = 1,
                                PERCENT_10_15 = 2,
                                PERCENT_15_20 = 3,
                                PERCENT_20_25 = 4,
                                PERCENT_25_30 = 5,
                                PERCENT_30_35 = 6,
                                PERCENT_35_40 = 7,
                                PERCENT_40_45 = 8,
                                PERCENT_45_50 = 9,
                                PERCENT_50_100 = 10};

constexpr static const size_t FWS_FREQUENCY_ARRAY_SIZE = 11;
using FwsFrequencyArray = std::array<AlleleSummmary, FWS_FREQUENCY_ARRAY_SIZE>;
using GenomeFWSMap = std::map<GenomeId_t, FwsFrequencyArray>;
using VariantFWSMap = std::map<std::string, AlleleSummmary>;

class CalcFWS {

public:

  CalcFWS() = default;
  ~CalcFWS() = default;

  void calcFwsStatistics(const std::shared_ptr<const PopulationDB>& population);
  [[nodiscard]] const GenomeFWSMap& getGenomeMap() const { return genome_fws_map_; }
  [[nodiscard]] const VariantFWSMap& getVariantMap() const { return variant_fws_map_; }
  void writeGenomeResults(const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr, const std::string& file_name) const;
  void writeVariantResults(const std::string& file_name) const;

private:

  GenomeFWSMap genome_fws_map_;
  VariantFWSMap variant_fws_map_;

  static std::pair<double, double> getFrequency(AlleleFrequencyBins bin_type);
  void updateGenomeFWSMap(const std::shared_ptr<const PopulationDB>& freq_population, size_t freq_bin);
  void updateVariantFWSMap(const std::shared_ptr<const PopulationDB>& population);

  constexpr static const char CSV_DELIMITER_ = ',';

};



} // Namespace.


#endif //KGL_ANALYSIS_PFEMP_FWS_H
