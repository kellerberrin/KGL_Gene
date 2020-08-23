//
// Created by kellerberrin on 23/8/20.
//

#include <kgl_variant_factory_vcf_evidence_analysis.h>
#include "kgl_variant.h"
#include "kgl_analysis_mutation_inbreed.h"


namespace kgl = kellerberrin::genome;


std::tuple<bool, double> kgl::InbreedingAnalysis::processFloatField( const std::shared_ptr<const Variant>& variant_ptr
                                                                   , const std::string& field_name) const {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, field_name);

  if (field_opt) {

    std::vector<double> field_vec = InfoEvidenceAnalysis::varianttoFloats(field_opt.value());

    if (field_vec.size() == 1) {

      return {true, field_vec.front() };

    } else if (field_vec.size() == 0) {

      // Missing value
      return {false, 0.0};

    } else {

      std::string vector_str;
      for (auto const& str : field_vec) {

        vector_str += str;
        vector_str += ";";

      }

      ExecEnv::log().error("InbreedingAnalysis::processField, Field: {} expected vector size 1, get vector size: {}, vector: {}",
                           field_name, field_vec.size(), vector_str);
      return {false, 0.0};

    }

  } else {

    ExecEnv::log().error("InbreedingAnalysis::processField, Field: {} not found for Variant: {}",
                         field_name, variant_ptr->output(',',VariantOutputIndex::START_0_BASED, false));

    return {false, 0.0};

  }

}



std::tuple<bool, double> kgl::InbreedingAnalysis::processStringField( const std::shared_ptr<const Variant>& variant_ptr,
                                                                      const std::string& field_name) const {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, field_name);

  if (field_opt) {

    std::vector<std::string> field_vec = InfoEvidenceAnalysis::varianttoStrings(field_opt.value());

    if (field_vec.size() != 1) {

      ExecEnv::log().error("InbreedingAnalysis::processStringField, Field: {} expected vector size 1, get vector size: {}",
                           field_name, field_vec.size());
      return {false, 0.0};

    } else {

      try {

        double frequency = std::stod(field_vec.front());
        return {true, frequency };

      }
      catch(std::exception& e) {

        ExecEnv::log().error("InbreedingAnalysis::processStringField, Field: {} problem converting: {} to double, exception: {}",
                             field_name, field_vec.front(), e.what());
        return {false, 0.0};

      }

    }

  } else {

    ExecEnv::log().error("InbreedingAnalysis::processField, Field: {} not found for Variant: {}",
                         field_name, variant_ptr->output(',',VariantOutputIndex::START_0_BASED, false));
    return {false, 0.0};
  }

}



std::string kgl::InbreedingAnalysis::lookupSuperPopulationField(const std::string& super_population) const {

  if (super_population == SUPER_POP_AFR_GNOMAD_.first) {

    return SUPER_POP_AFR_GNOMAD_.second;

  } else if (super_population == SUPER_POP_AMR_GNOMAD_.first) {

    return SUPER_POP_AMR_GNOMAD_.second;

  } else if (super_population == SUPER_POP_EAS_GNOMAD_.first) {

    return SUPER_POP_EAS_GNOMAD_.second;

  } else if (super_population == SUPER_POP_EUR_GNOMAD_.first) {

    return SUPER_POP_EUR_GNOMAD_.second;

  } else if (super_population == SUPER_POP_SAS_GNOMAD_.first) {

    return SUPER_POP_SAS_GNOMAD_.second;

  } else  {

    ExecEnv::log().error("MutationAnalysis::lookupSuperPopulationField; Unknown Super Population: {}", super_population);
    return SUPER_POP_SAS_GNOMAD_.second;

  }

}





// Join Diploid and a single genome population such as Gnomad or Clinvar.
void kgl::InbreedingAnalysis::joinPopulations() const {


  JoinSingleGenome join_pop(diploid_population_, unphased_population_, ped_data_);

  join_pop.joinPopulations();

  ExecEnv::log().info("Joined: {} Genomes in Population: {}", diploid_population_->getMap().size(), diploid_population_->populationId());

}



bool kgl::JoinSingleGenome::joinPopulations() {

  // The unphased population only has 1 genome. We find the variant by offset.
  if (joining_population_->getMap().size() != 1) {

    ExecEnv::log().error("JoinSingleGenome::lookupJoinedPop, unphased population : {} expected 1 genome, actually contains: {}",
                         joining_population_->populationId(), joining_population_->getMap().size());

    return false;

  }

  auto [joining_genome_id, joining_genome_ptr] = *(joining_population_->getMap().begin());

  ThreadPool thread_pool;
  std::vector<std::future<std::tuple<std::string, size_t, size_t>>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [joined_genome, joined_genome_ptr] : joined_population_->getMap()) {

    // function, object_ptr, arg1, ..., argn
    std::future<std::tuple<std::string, size_t, size_t>> future = thread_pool.enqueueTask(&JoinSingleGenome::processGenome, this, joined_genome_ptr, joining_genome_ptr);
    future_vector.push_back(std::move(future));

  }

  for (auto& future : future_vector) {

    auto [genome_id, total_variants, joined_variants] = future.get();

    double percent = (static_cast<double>(joined_variants) / static_cast<double>(total_variants)) * 100.0;

    std::string population_desc;
    auto result = ped_data_->getMap().find(genome_id);

    if (result == ped_data_->getMap().end()) {

      ExecEnv::log().warn("InbreedingAnalysis::checkPED, Genome sample: {} does not have a PED record", genome_id);
      population_desc = "";

    } else {

      population_desc = result->second.population();

    }

    ExecEnv::log().info("Population: {}, Genome: {}, Total Variants: {}, Joined Variants: {} ({}%)",
                        population_desc, genome_id, total_variants, joined_variants, percent);

  }

  return true;

}


std::tuple<std::string, size_t, size_t>
kgl::JoinSingleGenome::processGenome(std::shared_ptr<const DiploidGenome> diploid_genome, std::shared_ptr<const GenomeVariant> unphased_genome_ptr) {

  LocalGenomeJoin genome_join(unphased_genome_ptr);

  diploid_genome->processAll(genome_join, &LocalGenomeJoin::lookupJoinedPop);

  return {diploid_genome->genomeId(), genome_join.variantsProcessed(), genome_join.joinedVariantsFound()};

}

// Check that all samples (genomes) have a corresponding PED record.
void kgl::InbreedingAnalysis::checkPED() const {

  size_t PED_record_count = 0;
  for (auto const& [genome_id, genome_ptr] : diploid_population_->getMap()) {

    auto result = ped_data_->getMap().find(genome_id);

    if (result == ped_data_->getMap().end()) {

      ExecEnv::log().warn("InbreedingAnalysis::checkPED, Genome sample: {} does not have a PED record", genome_id);

    } else {

      ++PED_record_count;

    }

  }

  ExecEnv::log().info("Genome samples with PED records: {}", PED_record_count);

}


// Joins a single genome population (Gnomad, Clinvar) to another (generally phased Diploid) population.
bool kgl::LocalGenomeJoin::lookupJoinedPop(std::shared_ptr<const Variant> variant_ptr) {

  ++variants_processed_;

  auto contig_opt = unphased_genome_ptr_->getContig(variant_ptr->contigId());

  if (contig_opt) {

    auto variant_opt = contig_opt.value()->findVariant(*variant_ptr);

    if (variant_opt) {

      ++joined_variants_found_;

    }

  }

  return true;

}

