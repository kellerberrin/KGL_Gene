//
// Created by kellerberrin on 5/6/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H

#include "kgl_variant_factory_vcf_evidence.h"
#include "kgl_variant.h"


namespace kellerberrin::genome {   //  organization level namespace


// Utility Functions for extracting Info Data.
class InfoEvidenceAnalysis {

public:

  InfoEvidenceAnalysis() = default;
  ~InfoEvidenceAnalysis() = default;

  // The variant Info data (if it exists).
  static std::optional<const InfoSubscribedField> getSubscribedField( const std::shared_ptr<const Variant>& variant_ptr,
                                                                      const std::string& field_ident);

  static std::optional<InfoDataVariant> getInfoData( const std::shared_ptr<const Variant>& variant_ptr,
                                                     const std::string& field_ident);

  // Transform the returned data variant
  static std::vector<std::string> varianttoStrings(const InfoDataVariant& info_data);
  static std::vector<double> varianttoFloats(const InfoDataVariant& info_data);
  static std::vector<int64_t> varianttoIntegers(const InfoDataVariant& info_data);
  static bool variantToBool(const InfoDataVariant& info_data);

  // Converts a bin in string format "1|0|0|0|1|0|0|0|1|0" into a vector of floats.
  static std::vector<double> stringBinToFloat(const std::vector<std::string>& bin_data, size_t expected_bin_size);


public:

  constexpr static const char BIN_DELIMITER_{'|'};



};



} // namespace.


#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H
