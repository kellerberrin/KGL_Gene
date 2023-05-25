//
// Created by kellerberrin on 22/05/23.
//

#ifndef KGL_ANALYSIS_PFEMP_FWS_H
#define KGL_ANALYSIS_PFEMP_FWS_H


#include "kgl_analysis_PfEMP_heterozygous.h"


namespace kellerberrin::genome {   //  organization::project level namespace

enum class AlleleFrequencyBins: size_t { PERCENT_0_5 = 0,
                                PERCENT_5_10 = 1,
                                PERCENT_10_15 = 2,
                                PERCENT_15_20 = 3,
                                PERCENT_20_25 = 4,
                                PERCENT_25_30 = 5,
                                PERCENT_30_35 = 6,
                                PERCENT_35_40 = 7,
                                PERCENT_40_45 = 8,
                                PERCENT_45_50 = 9 };

constexpr static const size_t FWS_FREQUENCY_ARRAY_SIZE = 10;
using FwsFrequencyArray = std::array<VariantAnalysisType, FWS_FREQUENCY_ARRAY_SIZE>;
using GenomeFWSMap = std::map<GenomeId_t, FwsFrequencyArray>;

class CalcFWS {

public:

  CalcFWS() = default;
  ~CalcFWS() = default;

  void calcFwsStatistics(const std::shared_ptr<const PopulationDB>& population);
  [[nodiscard]] const GenomeFWSMap& getMap() const { return genome_fws_map_; }
  void writeResults(const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr, const std::string& file_name) const;

private:

  GenomeFWSMap genome_fws_map_;

  static std::pair<double, double> getFrequency(AlleleFrequencyBins bin_type);
  void updateGenomeFWSMap(const std::shared_ptr<const PopulationDB>& freq_population, size_t freq_bin);
  static void updateFreqRecord(const std::shared_ptr<const GenomeDB>& genome_ptr, VariantAnalysisType& freq_record);

  constexpr static const char CSV_DELIMITER_ = ',';

};



} // Namespace.


#endif //KGL_ANALYSIS_PFEMP_FWS_H
