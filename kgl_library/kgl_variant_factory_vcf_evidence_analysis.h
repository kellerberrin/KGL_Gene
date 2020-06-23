//
// Created by kellerberrin on 5/6/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H

#include "kgl_variant_factory_vcf_evidence.h"
#include "kgl_variant_db_unphased_population.h"


namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold parsed "vep" sub fields.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VEPSubFieldEvidence {

public:

  VEPSubFieldEvidence( std::shared_ptr<const VEPSubFieldHeader> vep_header_ptr,
                       std::vector<std::string>&& field_vector,
                       std::vector<std::vector<std::string_view>>&& vep_sub_fields_vector)
  : vep_header_ptr_(std::move(vep_header_ptr)),
    field_vector_(std::move(field_vector)),
    vep_sub_fields_vector_(std::move(vep_sub_fields_vector)) {}
  ~VEPSubFieldEvidence() = default;


  [[nodiscard]] std::shared_ptr<const VEPSubFieldHeader> vepHeader() const { return vep_header_ptr_; }
  [[nodiscard]] const std::vector<std::string>&  vepFields() const { return field_vector_; }
  [[nodiscard]] const std::vector<std::vector<std::string_view>>& vepSubFields() const { return vep_sub_fields_vector_; }

private:

  std::shared_ptr<const VEPSubFieldHeader> vep_header_ptr_;
  const std::vector<std::string> field_vector_;
  std::vector<std::vector<std::string_view>> vep_sub_fields_vector_;


};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Utility Functions for extracting Info Data.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoEvidenceAnalysis {

public:

  InfoEvidenceAnalysis() = default;
  ~InfoEvidenceAnalysis() = default;

  // The variant Info data (if it exists).
  static std::optional<const InfoSubscribedField> getSubscribedField( const Variant& variant,
                                                                      const std::string& field_ident);

  static std::optional<InfoDataVariant> getInfoData( const Variant& variant_ptr,
                                                     const std::string& field_ident);

  // Transform the returned data variant
  static const std::vector<std::string>& varianttoStrings(const InfoDataVariant& info_data);
  static const std::vector<double>& varianttoFloats(const InfoDataVariant& info_data);
  static const std::vector<int64_t>& varianttoIntegers(const InfoDataVariant& info_data);
  static bool variantToBool(const InfoDataVariant& info_data);

  // Converts a bin in string format "1|0|0|0|1|0|0|0|1|0" into a vector of floats.
  static std::vector<double> stringBinToFloat(const std::vector<std::string>& bin_data, size_t expected_bin_size);

  static std::optional<std::unique_ptr<const VEPSubFieldEvidence>> getVepSubFields(const Variant& variant);

  static void vepSubFieldValues( std::string vep_sub_field, const std::shared_ptr<const UnphasedPopulation>& population);


public:

  constexpr static const char BIN_DELIMITER_{'|'};



};






} // namespace.


#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H
