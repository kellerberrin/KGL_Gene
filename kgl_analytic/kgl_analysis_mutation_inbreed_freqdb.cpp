//
// Created by kellerberrin on 20/11/20.
//

#include <kgl_variant_factory_vcf_evidence_analysis.h>
#include "kgl_analysis_mutation_inbreed_freqdb.h"


namespace kgl = kellerberrin::genome;

std::optional<double> kgl::FrequencyDatabaseRead::processFloatField(const Variant& variant,
                                                                    const std::string& frequency_field) const {

  // Use a super population code to lookup a corresponding AF field.
  std::string database_field = lookupVariantSuperPopField(frequency_field);

  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(variant, database_field);

  if (field_opt) {

    std::vector<double> field_vec = InfoEvidenceAnalysis::varianttoFloats(field_opt.value());

    if (field_vec.size() == 1) {

      return field_vec.front();

    } else if (field_vec.size() == 0) {

      // Missing value, this is OK, means the field exists but not defined for this variant.
      return std::nullopt;

    } else {

      std::string vector_str;
      for (auto const& str : field_vec) {

        vector_str += str;
        vector_str += ";";

      }

      ExecEnv::log().warn("InbreedingAnalysis::processFloatField, Field: {} expected vector size 1, get vector size: {}, vector: {}, Variant: {}",
                          database_field, field_vec.size(), vector_str, variant.output(',', VariantOutputIndex::START_0_BASED, false));
      return std::nullopt;

    }

  } else {

    ExecEnv::log().warn("InbreedingAnalysis::processFloatField, Field: {} Not found for variant: {}",
                        database_field, variant.output(',', VariantOutputIndex::START_0_BASED, false));
    return std::nullopt;

  }

}


std::string kgl::FrequencyDatabaseRead::lookupVariantSuperPopField(const std::string& super_population) const {

  switch(source_) {

    case FrequencyDatabaseSource::GNOMAD2_1:
      return lookupGnomad_2_1_Field(super_population);

    case FrequencyDatabaseSource::GNOMAD3_1:
      return lookupGnomad_3_1_Field(super_population);

    case FrequencyDatabaseSource::GENOMES_1000:
      return lookup_1000_Field(super_population);

  }

  // To placate the compiler.
  return lookupGnomad_2_1_Field(super_population);

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

    ExecEnv::log().error("MutationAnalysis::lookupGnomad_2_1_Field; Unknown Super Population: {}", super_population);
    return SUPER_POP_ALL_GNOMAD_2_1.second;

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

    ExecEnv::log().error("MutationAnalysis::lookupGnomad_3_1_Field; Unknown Super Population: {}", super_population);
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

    ExecEnv::log().error("MutationAnalysis::lookup_1000_Field; Unknown Super Population: {}", super_population);
    return SUPER_POP_SAS_1000.second;

  }

}
