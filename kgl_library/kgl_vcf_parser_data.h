//
// Created by kellerberrin on 8/04/18.
//

#ifndef KGL_VCF_PARSER_DATA_H
#define KGL_VCF_PARSER_DATA_H


#include "kgl_variant_db.h"
#include "kgl_variant_phasing_statistics.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// General  Parser Analysis
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class ParserAnalysis : public PopulationVariant {


public:

  ParserAnalysis(const std::string& analysis) : PopulationVariant(analysis),
                                                phased_statistics_ptr(std::make_shared<PopulationPhasingStatistics>()) {}
  ~ParserAnalysis() override = default;

  std::shared_ptr<PopulationPhasingStatistics> phasedStatistics() { return phased_statistics_ptr; }
  std::shared_ptr<PopulationPhasingStatistics> phasedStatistics() const { return phased_statistics_ptr; }


private:

  std::shared_ptr<PopulationPhasingStatistics> phased_statistics_ptr;

};






}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_KGL_VCF_PARSER_DATA_H
