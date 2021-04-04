//
// Created by kellerberrin on 20/11/20.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_variant_db_freq.h"


namespace kgl = kellerberrin::genome;



std::optional<double> kgl::FrequencyDatabaseRead::processSuperPopField(const Variant& variant, const std::string& super_population) {

  // Use a super population code to lookup a corresponding AF field.
  std::optional<std::string> database_field_opt = lookupVariantSuperPopField(variant.evidence().dataSource(), super_population);

  if (not database_field_opt) {

    ExecEnv::log().warn("FrequencyDatabaseRead::processSuperPopField; Unable to to lookup superpopulation frequency for data source");
    return 0.0;

  }

  auto database_field = database_field_opt.value();

  return processFloatField(variant, database_field);

}


std::optional<double> kgl::FrequencyDatabaseRead::processFloatField(const Variant& variant, const std::string& database_field) {

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

      ExecEnv::log().warn("FrequencyDatabaseRead::processFloatField; Field: {} expected vector size 1, evidence variants: {}, evidence index: {},  get vector size: {}, vector: {}, Variant: {}",
                          database_field,
                          variant.evidence().altVariantCount(),
                          variant.evidence().altVariantIndex(),
                          field_vec.size(),
                          vector_str,
                          variant.output(',', VariantOutputIndex::START_0_BASED, true));
      return std::nullopt;

    }

  } else {

    ExecEnv::log().warn("FrequencyDatabaseRead::processFloatField; Field: {} Not found for variant: {}",
                        database_field, variant.output(',', VariantOutputIndex::START_0_BASED, true));
    return std::nullopt;

  }

}





std::optional<std::string> kgl::FrequencyDatabaseRead::lookupVariantSuperPopField( DataSourceEnum data_source,
                                                                                   const std::string& super_population) {

  switch(data_source) {

    case DataSourceEnum::Gnomad2_1:
      return lookupGnomad_2_1_Field(super_population);

    case DataSourceEnum::GnomadExomes2_1:
      return lookupGnomad_ex_2_1_Field(super_population);

    case DataSourceEnum::Gnomad3_1:
      return lookupGnomad_3_1_Field(super_population);

    case DataSourceEnum::GnomadGenome3_1:
      return lookupGnomadGenome_3_1_Field(super_population);

    case DataSourceEnum::Genome1000:
      return lookup_1000_Field(super_population);

    default:
      return std::nullopt;

  }

}



std::string kgl::FrequencyDatabaseRead::lookupGnomad_2_1_Field(const std::string& super_population) {

  if (super_population == SUPER_POP_AFR_GNOMAD_2_1.first) {

    return SUPER_POP_AFR_GNOMAD_2_1.second;

  } else if (super_population == SUPER_POP_AMR_GNOMAD_2_1.first) {

    return SUPER_POP_AMR_GNOMAD_2_1.second;

  } else if (super_population == SUPER_POP_EAS_GNOMAD_2_1.first) {

    return SUPER_POP_EAS_GNOMAD_2_1.second;

  } else if (super_population == SUPER_POP_EUR_GNOMAD_2_1.first) {

    return SUPER_POP_EUR_GNOMAD_2_1.second;

  } else if (super_population == SUPER_POP_SAS_GNOMAD_2_1.first) {

    return SUPER_POP_SAS_GNOMAD_2_1.second;

  } else if (super_population == SUPER_POP_ALL_GNOMAD_2_1.first) {

    return SUPER_POP_ALL_GNOMAD_2_1.second;

  } else  {

    ExecEnv::log().error("FrequencyDatabaseRead::lookupGnomad_2_1_Field; Unknown Super Population: {}", super_population);
    return SUPER_POP_ALL_GNOMAD_2_1.second;

  }

}


std::string kgl::FrequencyDatabaseRead::lookupGnomad_ex_2_1_Field(const std::string& super_population) {

  if (super_population == SUPER_POP_AFR_GNOMAD_2_1.first) {

    return SUPER_POP_AFR_GNOMAD_EX_2_1.second;

  } else if (super_population == SUPER_POP_AMR_GNOMAD_2_1.first) {

    return SUPER_POP_AMR_GNOMAD_EX_2_1.second;

  } else if (super_population == SUPER_POP_EAS_GNOMAD_2_1.first) {

    return SUPER_POP_EAS_GNOMAD_EX_2_1.second;

  } else if (super_population == SUPER_POP_EUR_GNOMAD_2_1.first) {

    return SUPER_POP_EUR_GNOMAD_EX_2_1.second;

  } else if (super_population == SUPER_POP_SAS_GNOMAD_2_1.first) {

    return SUPER_POP_SAS_GNOMAD_EX_2_1.second;

  } else if (super_population == SUPER_POP_ALL_GNOMAD_2_1.first) {

    return SUPER_POP_ALL_GNOMAD_EX_2_1.second;

  } else  {

    ExecEnv::log().error("FrequencyDatabaseRead::lookupGnomad_ex-2_1_Field; Unknown Super Population: {}", super_population);
    return SUPER_POP_ALL_GNOMAD_EX_2_1.second;

  }

}


std::string kgl::FrequencyDatabaseRead::lookupGnomad_3_1_Field(const std::string& super_population) {

  if (super_population == SUPER_POP_AFR_GNOMAD_3_1.first) {

    return SUPER_POP_AFR_GNOMAD_3_1.second;

  } else if (super_population == SUPER_POP_AMR_GNOMAD_3_1.first) {

    return SUPER_POP_AMR_GNOMAD_3_1.second;

  } else if (super_population == SUPER_POP_EAS_GNOMAD_3_1.first) {

    return SUPER_POP_EAS_GNOMAD_3_1.second;

  } else if (super_population == SUPER_POP_EUR_GNOMAD_3_1.first) {

    return SUPER_POP_EUR_GNOMAD_3_1.second;

  } else if (super_population == SUPER_POP_SAS_GNOMAD_3_1.first) {

    return SUPER_POP_SAS_GNOMAD_3_1.second;

  } else if (super_population == SUPER_POP_ALL_GNOMAD_3_1.first) {

    return SUPER_POP_ALL_GNOMAD_3_1.second;

  } else  {

    ExecEnv::log().error("FrequencyDatabaseRead::lookupGnomad_3_1_Field; Unknown Super Population: {}", super_population);
    return SUPER_POP_ALL_GNOMAD_3_1.second;

  }

}

std::string kgl::FrequencyDatabaseRead::lookupGnomadGenome_3_1_Field(const std::string& super_population) {

  if (super_population == SUPER_POP_AFR_GNOMAD_3_1.first) {

    return SUPER_POP_AFR_GNOMADGENOME_3_1.second;

  } else if (super_population == SUPER_POP_AMR_GNOMAD_3_1.first) {

    return SUPER_POP_AMR_GNOMADGENOME_3_1.second;

  } else if (super_population == SUPER_POP_EAS_GNOMAD_3_1.first) {

    return SUPER_POP_EAS_GNOMADGENOME_3_1.second;

  } else if (super_population == SUPER_POP_EUR_GNOMAD_3_1.first) {

    return SUPER_POP_EUR_GNOMADGENOME_3_1.second;

  } else if (super_population == SUPER_POP_SAS_GNOMAD_3_1.first) {

    return SUPER_POP_SAS_GNOMADGENOME_3_1.second;

  } else if (super_population == SUPER_POP_ALL_GNOMAD_3_1.first) {

    return SUPER_POP_ALL_GNOMADGENOME_3_1.second;

  } else  {

    ExecEnv::log().error("FrequencyDatabaseRead::lookupGnomadGenome_3_1_Field; Unknown Super Population: {}", super_population);
    return SUPER_POP_ALL_GNOMAD_3_1.second;

  }

}



std::string kgl::FrequencyDatabaseRead::lookup_1000_Field(const std::string& super_population) {

  if (super_population == SUPER_POP_AFR_1000.first) {

    return SUPER_POP_AFR_1000.second;

  } else if (super_population == SUPER_POP_AMR_1000.first) {

    return SUPER_POP_AMR_1000.second;

  } else if (super_population == SUPER_POP_EAS_1000.first) {

    return SUPER_POP_EAS_1000.second;

  } else if (super_population == SUPER_POP_EUR_1000.first) {

    return SUPER_POP_EUR_1000.second;

  } else if (super_population == SUPER_POP_SAS_1000.first) {

    return SUPER_POP_SAS_1000.second;

  } else if (super_population == SUPER_POP_ALL_1000.first) {

    return SUPER_POP_ALL_1000.second;

  } else  {

    ExecEnv::log().error("FrequencyDatabaseRead::lookup_1000_Field; Unknown Super Population: {}", super_population);
    return SUPER_POP_SAS_1000.second;

  }

}