//
// Created by kellerberrin on 20/11/20.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_variant_db_freq.h"


namespace kgl = kellerberrin::genome;



// Map a data source to a super population field code index. Returns false if the data source is unsupported.
bool kgl::FrequencyDatabaseRead::sourceIndex(DataSourceEnum data_source, size_t& source_index) {

  switch (data_source) {

    case DataSourceEnum::Gnomad2_1:
      source_index = 0;
      return true;

    case DataSourceEnum::GnomadExomes2_1:
      source_index = 1;
      return true;

    case DataSourceEnum::Gnomad3_1:
      source_index = 2;
      return true;

    case DataSourceEnum::GnomadGenome3_1:
      source_index = 3;
      return true;

    case DataSourceEnum::Genome1000:
      source_index = 4;
      return true;

    default:
      return false;

  }

}


// Lookup a super population field code. Returns nullopt if the super population or data source is unsupported.
std::optional<std::string_view> kgl::FrequencyDatabaseRead::lookupSuperPopField( FieldIndex field_index,
                                                                                 DataSourceEnum data_source,
                                                                                 const std::string& super_population) {

  size_t source_index{0};
  if (not sourceIndex(data_source, source_index)) {

    ExecEnv::log().warn("FrequencyDatabaseRead; Unable to lookup superpopulation field, unsupported data source: {}",
                        DataDB::dataSource(data_source));
    return std::nullopt;

  }

  // Select the AF, AN or AC field code table.
  const SuperPopFields& field_table = field_text_tables_[static_cast<size_t>(field_index)];

  auto field_iter = std::ranges::find(field_table, super_population, [](const SuperPopCode& field_pair) { return field_pair.first; });
  if (field_iter == field_table.end()) {

    ExecEnv::log().warn("FrequencyDatabaseRead; Unsupported super population: {}", super_population);
    return std::nullopt;

  }

  return field_iter->second[source_index];

}


std::optional<double> kgl::FrequencyDatabaseRead::superPopFrequency(const Variant& variant, const std::string& super_population) {

  // Use a super population code to lookup a corresponding AF field.
  auto database_field_opt = lookupSuperPopField(FieldIndex::AF, variant.evidence().dataSource(), super_population);

  if (not database_field_opt) {

    return std::nullopt;

  }

  return infoFloatField(variant, std::string{database_field_opt.value()});

}


std::optional<int64_t> kgl::FrequencyDatabaseRead::superPopTotalAlleles(const Variant& variant, const std::string& super_population) {

  // Use a super population code to lookup a corresponding AN field.
  auto database_field_opt = lookupSuperPopField(FieldIndex::AN, variant.evidence().dataSource(), super_population);

  if (not database_field_opt) {

    return std::nullopt;

  }

  return infoIntegerField(variant, std::string{database_field_opt.value()});

}


std::optional<int64_t> kgl::FrequencyDatabaseRead::superPopAltAlleles(const Variant& variant, const std::string& super_population) {

  // Use a super population code to lookup a corresponding AC field.
  auto database_field_opt = lookupSuperPopField(FieldIndex::AC, variant.evidence().dataSource(), super_population);

  if (not database_field_opt) {

    return std::nullopt;

  }

  return infoIntegerField(variant, std::string{database_field_opt.value()});

}


// Common implementation for infoFloatField() and infoIntegerField().
template<typename T>
std::optional<T> kgl::FrequencyDatabaseRead::infoFieldImpl( const Variant& variant,
                                                            const std::string& database_field,
                                                            const char* field_type) {

  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(variant, database_field);

  if (not field_opt) {

    ExecEnv::log().warn("FrequencyDatabaseRead::info{}Field; Field: {} Not found for variant: {}",
                        field_type, database_field, variant.HGVS_Phase());
    return std::nullopt;

  }

  // Convert the info data to a vector of the requested type (double or int64_t).
  const std::vector<T>& field_vec = [&]() -> const std::vector<T>& {

    if constexpr (std::is_same_v<T, double>) {

      return InfoEvidenceAnalysis::varianttoFloats(field_opt.value());

    } else {

      return InfoEvidenceAnalysis::varianttoIntegers(field_opt.value());

    }

  }();

  if (field_vec.empty()) {

    // Missing value, this is OK, means the field exists but not defined for this variant.
    return std::nullopt;

  } else if (field_vec.size() == 1) {

    return field_vec.front();

  } else if (variant.evidence().altVariantCount() == field_vec.size()
             and variant.evidence().altVariantIndex() < field_vec.size()) {

    return field_vec[variant.evidence().altVariantIndex()];

  }

  // Vector size mismatch, log and return missing.
  std::string vector_str;
  for (auto const& value : field_vec) {

    vector_str += std::to_string(value);
    vector_str += ";";

  }

  ExecEnv::log().warn("FrequencyDatabaseRead::info{}Field; Field: {} expected vector size 1, evidence variants: {}, evidence index: {},  get vector size: {}, vector: {}, Variant: {}",
                      field_type,
                      database_field,
                      variant.evidence().altVariantCount(),
                      variant.evidence().altVariantIndex(),
                      field_vec.size(),
                      vector_str,
                      variant.HGVS_Phase());
  return std::nullopt;

}


std::optional<double> kgl::FrequencyDatabaseRead::infoFloatField(const Variant& variant, const std::string& database_field) {

  return infoFieldImpl<double>(variant, database_field, "Float");

}


std::optional<int64_t> kgl::FrequencyDatabaseRead::infoIntegerField(const Variant& variant, const std::string& database_field) {

  return infoFieldImpl<int64_t>(variant, database_field, "Integer");

}