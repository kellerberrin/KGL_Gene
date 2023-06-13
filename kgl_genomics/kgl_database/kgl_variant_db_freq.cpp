//
// Created by kellerberrin on 20/11/20.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_variant_db_freq.h"


namespace kgl = kellerberrin::genome;



std::optional<double> kgl::FrequencyDatabaseRead::superPopFrequency(const Variant& variant, const std::string& super_population) {

  // Use a super population code to lookup a corresponding AF field.
  std::optional<std::string> database_field_opt = lookupVariantSuperPopField(field_text_map_AF_, variant.evidence().dataSource(), super_population);

  if (not database_field_opt) {

    ExecEnv::log().warn("FrequencyDatabaseRead::superPopFrequency; Unable to to lookup superpopulation frequency for data source");
    return 0.0;

  }

  auto database_field = database_field_opt.value();

  return infoFloatField(variant, database_field);

}



std::optional<int64_t> kgl::FrequencyDatabaseRead::superPopTotalAlleles(const Variant& variant, const std::string& super_population) {

  // Use a super population code to lookup a corresponding AF field.
  std::optional<std::string> database_field_opt = lookupVariantSuperPopField(field_text_map_AN_, variant.evidence().dataSource(), super_population);

  if (not database_field_opt) {

    ExecEnv::log().warn("FrequencyDatabaseRead::superPopTotalAlleles; Unable to to lookup superpopulation frequency for data source");
    return 0.0;

  }

  auto database_field = database_field_opt.value();

  return infoIntegerField(variant, database_field);

}


std::optional<int64_t> kgl::FrequencyDatabaseRead::superPopAltAlleles(const Variant& variant, const std::string& super_population) {

  // Use a super population code to lookup a corresponding AF field.
  std::optional<std::string> database_field_opt = lookupVariantSuperPopField(field_text_map_AC_, variant.evidence().dataSource(), super_population);

  if (not database_field_opt) {

    ExecEnv::log().warn("FrequencyDatabaseRead::superPopTotalAlleles; Unable to to lookup superpopulation frequency for data source");
    return 0.0;

  }

  auto database_field = database_field_opt.value();

  return infoIntegerField(variant, database_field);

}



std::optional<double> kgl::FrequencyDatabaseRead::infoFloatField(const Variant& variant, const std::string& database_field) {

  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(variant, database_field);

  if (field_opt) {

    std::vector<double> field_vec = InfoEvidenceAnalysis::varianttoFloats(field_opt.value());

    if (field_vec.size() == 1) {

      return field_vec.front();

    } else if (field_vec.empty()) {

      // Missing value, this is OK, means the field exists but not defined for this variant.
      return std::nullopt;

    } else if (variant.evidence().altVariantCount() == field_vec.size()
           and variant.evidence().altVariantIndex() < field_vec.size()) {

      return field_vec[variant.evidence().altVariantIndex()];

    } else {

      std::string vector_str;
      for (auto const& str : field_vec) {

        vector_str += std::to_string(str);
        vector_str += ";";

      }

      ExecEnv::log().warn("FrequencyDatabaseRead::infoFloatField; Field: {} expected vector size 1, evidence variants: {}, evidence index: {},  get vector size: {}, vector: {}, Variant: {}",
                          database_field,
                          variant.evidence().altVariantCount(),
                          variant.evidence().altVariantIndex(),
                          field_vec.size(),
                          vector_str,
                          variant.HGVS_Phase());
      return std::nullopt;

    }

  } else {

    ExecEnv::log().warn("FrequencyDatabaseRead::infoFloatField; Field: {} Not found for variant: {}", database_field, variant.HGVS_Phase());
    return std::nullopt;

  }

}



std::optional<int64_t> kgl::FrequencyDatabaseRead::infoIntegerField(const Variant& variant, const std::string& database_field) {

  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(variant, database_field);

  if (field_opt) {

    std::vector<int64_t> field_vec = InfoEvidenceAnalysis::varianttoIntegers(field_opt.value());

    if (field_vec.size() == 1) {

      return field_vec.front();

    } else if (field_vec.empty()) {

      // Missing value, this is OK, means the field exists but not defined for this variant.
      return std::nullopt;

    } else if (variant.evidence().altVariantCount() == field_vec.size()
               and variant.evidence().altVariantIndex() < field_vec.size()) {

      return field_vec[variant.evidence().altVariantIndex()];

    } else {

      std::string vector_str;
      for (auto const& str : field_vec) {

        vector_str += std::to_string(str);
        vector_str += ";";

      }

      ExecEnv::log().warn("FrequencyDatabaseRead::infoIntegerField; Field: {} expected vector size 1, evidence variants: {}, evidence index: {},  get vector size: {}, vector: {}, Variant: {}",
                          database_field,
                          variant.evidence().altVariantCount(),
                          variant.evidence().altVariantIndex(),
                          field_vec.size(),
                          vector_str,
                          variant.HGVS_Phase());
      return std::nullopt;

    }

  } else {

    ExecEnv::log().warn("FrequencyDatabaseRead::infoIntegerField; Field: {} Not found for variant: {}",
                        database_field, variant.HGVS_Phase());
    return std::nullopt;

  }

}


std::optional<std::string> kgl::FrequencyDatabaseRead::lookupVariantSuperPopField( const std::map<std::string, SuperPopFieldText>& lookup_map,
                                                                                   DataSourceEnum data_source,
                                                                                   const std::string& super_population) {

  auto result = lookup_map.find(super_population);
  if (result == lookup_map.end()) {

    return std::nullopt;

  }

  auto const& [super_pop, field_text] = *result;

  switch(data_source) {

    case DataSourceEnum::Gnomad2_1:
      return field_text.gnomad_2_1;

    case DataSourceEnum::GnomadExomes2_1:
      return field_text.gnomad_ex_2_1;

    case DataSourceEnum::Gnomad3_1:
      return field_text.gnomad_3_1;

    case DataSourceEnum::GnomadGenome3_1:
      return field_text.gnomadgenome_3_1;

    case DataSourceEnum::Genome1000:
      return field_text.genome_1000;

    default:
      return std::nullopt;

  }

}


