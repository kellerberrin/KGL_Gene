//
// Created by kellerberrin on 27/08/18.
//

#ifndef KGL_READ_PHASING_H
#define KGL_READ_PHASING_H

#include "kgl_genome_types.h"
#include <map>
#include <array>



namespace kellerberrin::genome {   //  organization::project level namespace


using MixtureStatistics = std::pair<size_t, std::array<double, 5>>;  // The sample mixture data. size_t is the strain count.
using MixtureStatisticsMap = std::map<GenomeId_t, MixtureStatistics>; // THe mixture stats indexed by genome.
class GenomeMixtureStatistics {

public:

  GenomeMixtureStatistics() = default;
  ~GenomeMixtureStatistics() = default;

  const MixtureStatisticsMap& getMap() const { return mixture_statistics_map_; }
  // Returns false if no genome found, mixture_statistics is zeroed.
  bool getMixtureStatistics(const GenomeId_t& genome, MixtureStatistics& mixture_statistics) const;
  // Read the mixture file. Can be more than 1 char used as field delimiter.
  bool readMixtureStatistics(const std::string& file_name, const std::string& field_delimiters = FIELD_DELIMITER_);

private:

  MixtureStatisticsMap mixture_statistics_map_;

  constexpr static const char* FIELD_DELIMITER_ = "\t";
  constexpr static const char COMMENT_CHAR_ = '#';
  constexpr static const size_t MIN_FIELD_COUNT_ = 7;
  constexpr static const double TOLERANCE_ = 0.0001;
  constexpr static const double SINGLE_STRAIN_ = 0.99;

};



}   // end namespace



#endif //KGL_READ_PHASING_H
