//
// Created by kellerberrin on 5/6/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H

#include "kgl_variant_factory_vcf_evidence.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold parsed "vep" sub fields.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VEPSubFieldEvidence {

public:

  VEPSubFieldEvidence( std::shared_ptr<const VEPSubFieldHeader> vep_header_ptr,
                       std::vector<std::string>&& vep_fields_vector) : vep_header_ptr_(std::move(vep_header_ptr)),
                                                                       vep_fields_vector_(std::move(vep_fields_vector)) {}
  ~VEPSubFieldEvidence() = default;

  [[nodiscard]] std::shared_ptr<const VEPSubFieldHeader> vepHeader() const { return vep_header_ptr_; }
  [[nodiscard]] const std::vector<std::string>& vepFields() const { return vep_fields_vector_; }
  // Adapter converts a vep field to sub_fields.
  [[nodiscard]] static const std::vector<std::string_view> vepSubFields(const std::string& vep_field);

private:

  std::shared_ptr<const VEPSubFieldHeader> vep_header_ptr_;
  const std::vector<std::string> vep_fields_vector_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Utility Functions for extracting Info Data.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Indexed by vep field.
using VepValueMap = std::vector<std::map<std::string, std::string>>;
using VepIndexVector = std::vector<std::pair<std::string, size_t>>;


template<typename T>
concept ValidInfoDataType =
std::same_as<T, std::vector<int64_t>> ||
std::same_as<T, std::vector<double>> ||
std::same_as<T, std::vector<std::string>> ||
                std::same_as<T, int64_t> ||
                std::same_as<T, double> ||
                std::same_as<T, std::string> ||
                std::same_as<T, bool>;

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

  // There can multiple vep records per variant. Default empty vector retrieves all fields.
  static VepValueMap getVepValues(const Variant& variant, std::vector<std::string> vep_fields = std::vector<std::string>{});
  // Get the vep field offsets. Typically this only need to be done once if all variants have the same vep layout.
  static VepIndexVector getVepIndexes(const Variant& variant, const std::vector<std::string>& vep_field_list);
  // Use the indexes generated above retrieve the vep data.
  static VepValueMap getVepData(const Variant& variant, const VepIndexVector& vep_field_list);

  // Display all discrete values.
  static void vepSubFieldValues( std::string vep_sub_field, const std::shared_ptr<const PopulationDB>& population);
  // Count discrete values.
  static void vepSubFieldCount( std::string vep_sub_field, const std::shared_ptr<const PopulationDB>& population);
  // Type template getTypedData
  template<typename T> requires ValidInfoDataType<T>
  static std::optional<T> getTypedInfoData( const Variant& variant, const std::string& field_ident);

public:

  constexpr static const char BIN_DELIMITER_{'|'};
  inline static size_t vep_size_errors_{0};
  inline static size_t vep_records_{0};



};

// Can only be templated with double, std::vector<double>, int64_t, std::vector<int64_t>, std::string,
// std::vector<string> and bool (add concepts).
template<typename T> requires ValidInfoDataType<T>
std::optional<T> InfoEvidenceAnalysis::getTypedInfoData( const Variant& variant, const std::string& field_ident) {

  std::optional<const InfoSubscribedField> field_opt = getSubscribedField(variant, field_ident);

  if (field_opt) {

    auto info_data_ptr = variant.evidence().infoData();
    if (info_data_ptr) {

      const DataMemoryBlock &data_block = *(info_data_ptr.value());

      {

        auto const& field_info = field_opt.value();
        // Check that this is the correct index for the data block.
        // Should never happen (but best to check).
        if (field_info.getDataHeader() != data_block.evidenceHeader()) {

          ExecEnv::log().error( "InfoEvidenceAnalysis::getTypedInfoData, Incorrect Subscribed Field Index used Access Info Data Block");
          return std::nullopt; // Inaccessible variant data, will trigger downstream errors.

        }

        auto resource_type = field_info.getDataHandle().resourceType();

        if constexpr(std::is_same_v<T, double>) {

          if (resource_type != DataResourceType::Float) {

            return std::nullopt;

          }

          auto info_data = data_block.getFloat(field_info.getDataHandle());
          if (info_data.size() != 1) {

            return std::nullopt;

          }

          return info_data.front();

        }

        if constexpr(std::is_same_v<T, std::vector<double>>) {

          if (resource_type != DataResourceType::Float) {

            return std::nullopt;

          }

          auto info_data = data_block.getFloat(field_info.getDataHandle());

          return info_data;

        }

        if constexpr(std::is_same_v<T, int64_t>) {

          if (resource_type != DataResourceType::Integer) {

            return std::nullopt;

          }

          auto info_data = data_block.getInteger(field_info.getDataHandle());
          if (info_data.size() != 1) {

            return std::nullopt;

          }

          return info_data.front();

        }

        if constexpr(std::is_same_v<T, std::vector<int64_t>>) {

          if (resource_type != DataResourceType::Integer) {

            return std::nullopt;

          }

          auto info_data = data_block.getInteger(field_info.getDataHandle());

          return info_data;

        }


        if constexpr(std::is_same_v<T, std::string>) {

          if (resource_type != DataResourceType::String) {

            return std::nullopt;

          }

          auto info_data = data_block.getString(field_info.getDataHandle());
          if (info_data.size() != 1) {

            return std::nullopt;

          }

          return info_data.front();

        }

        if constexpr(std::is_same_v<T, std::vector<std::string>>) {

          if (resource_type != DataResourceType::String) {

            return std::nullopt;

          }

          auto info_data = data_block.getInteger(field_info.getDataHandle());

          return info_data;

        }

        if constexpr(std::is_same_v<T, bool>) {

          if (resource_type != DataResourceType::Boolean) {

            return std::nullopt;

          }

          auto info_data = data_block.getBoolean(field_info.getDataHandle());

          return info_data;

        }

        ExecEnv::log().error( "InfoEvidenceAnalysis::getTypedInfoData, Internal data type unknown");

      }

    }

  }

  return std::nullopt;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using VepFieldValueMap = std::map<std::string, size_t>;

class VepSubFieldValues {

public:

  VepSubFieldValues(std::string vep_sub_field) : vep_sub_field_(std::move(vep_sub_field)) {}
  ~VepSubFieldValues() = default;

  [[nodiscard]] const VepFieldValueMap& getMap() const { return field_value_map_; }

  // For a single variant
  bool getSubFieldValues(const std::shared_ptr<const Variant>& variant_ptr);
  // For a population.
  bool getPopulationValues(const std::shared_ptr<const PopulationDB>& population_ptr);
  // For a genome
  bool getGenomeValues(const std::shared_ptr<const GenomeDB>& genome_ptr);
  // For a contig.
  bool getContigValues(const std::shared_ptr<const ContigDB>& contig_ptr);

private:

  const std::string vep_sub_field_;
  VepFieldValueMap field_value_map_;


}; // struct.




} // namespace.


#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_ANALYSIS_H
