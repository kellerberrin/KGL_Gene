//
// Created by kellerberrin on 5/6/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H

#include "kgl_variant_factory_vcf_evidence.h"


namespace kellerberrin::genome {   //  organization level namespace



class InfoEvidenceAnalysis {

  InfoEvidenceAnalysis() = default;
  ~InfoEvidenceAnalysis() = default;

  static std::vector<std::string> varianttoStrings(const InfoDataVariant&& info_data);
  static std::vector<double> varianttoFloats(const InfoDataVariant&& info_data);
  static std::vector<int64_t> varianttoIntegers(const InfoDataVariant&& info_data);
  static bool variantToBool(const InfoDataVariant&& info_data);

  static std::vector<double> stringBinToVector(const std::vector<std::string>& bin_data, size_t expected_bin_size = 10);


public:

  constexpr static const char BIN_DELIMITER_{'|'};



};


} // namespace.


#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H
