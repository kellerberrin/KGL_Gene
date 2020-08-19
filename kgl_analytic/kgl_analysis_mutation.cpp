//
// Created by kellerberrin on 3/7/20.
//

#include "kel_thread_pool.h"
#include "kgl_analysis_mutation.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"


#include <fstream>

namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::MutationAnalysis::initializeAnalysis(const std::string& work_directory,
                                               const RuntimeParameterMap& named_parameters,
                                               std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Default Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  std::optional<std::shared_ptr<const GenomeReference>> ref_genome_opt = reference_genomes->getOptionalGenome(REFERENCE_GENOME_);

  if (ref_genome_opt) {

    genome_GRCh38_ = ref_genome_opt.value();

  } else {

    ExecEnv::log().error("MutationAnalysis::initializeAnalysis, Could not find Genome: {} Analysis: {} disabled.", REFERENCE_GENOME_, ident());
    return false;

  }

  if (not getParameters(work_directory, named_parameters)) {

    return false;

  }

// Clear the data file a
  std::ofstream outfile;
  outfile.open(output_file_name_, std::ofstream::out | std::ofstream::trunc);

  return true;

}

// This function superclasses the data objects and stores them for further use.
bool kgl::MutationAnalysis::fileReadAnalysis(std::shared_ptr<const DataObjectBase> data_object_ptr) {

  ExecEnv::log().info("Analysis: {}, begin processing data file", ident(), data_object_ptr->Id());


  if (data_object_ptr->dataType() == DataTypeEnum::DiploidPopulation) {

    std::shared_ptr<const DiploidPopulation> diploid_population = std::dynamic_pointer_cast<const DiploidPopulation>(data_object_ptr);

    if (diploid_population) {

      ExecEnv::log().info("Analysis: {}, Generate Hetreozygous/Homozygous ratio statistics for file: {}", ident(), data_object_ptr->Id());

      diploid_population_ = diploid_population;

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a Diploid Population", ident(), data_object_ptr->Id());
      return false;

    }

  }

  if (data_object_ptr->dataType() == DataTypeEnum::UnphasedPopulation) {

    std::shared_ptr<const UnphasedPopulation> unphased_population = std::dynamic_pointer_cast<const UnphasedPopulation>(data_object_ptr);

    if (unphased_population) {

      unphased_population_ = unphased_population;

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not an Unphased Population", ident(), data_object_ptr->Id());
      return false;

    }

  }

  if (data_object_ptr->dataType() == DataTypeEnum::PedAncestor) {

    std::shared_ptr<const GenomePEDData> ped_data = std::dynamic_pointer_cast<const GenomePEDData>(data_object_ptr);

    if (ped_data) {

      ped_data_ = ped_data;
      ExecEnv::log().info("Analysis: {}, ped file: {} contains: {} PED records", ident(), ped_data->Id(), ped_data->getMap().size());

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a PED Ancestor Object", ident(), data_object_ptr->Id());
      return false;

    }

  }

  ExecEnv::log().info("Analysis: {}, completed data file: {}", ident(), data_object_ptr->Id());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Iteration Analysis called for Analysis Id: {}", ident());

  if (diploid_population_ and unphased_population_ and ped_data_) {

    ExecEnv::log().info("Filtered Population : {}  and Joined Population: {}  both active",
                        diploid_population_->populationId(), unphased_population_->populationId());

    if (not hetHomRatioLocus(diploid_population_)) {

      ExecEnv::log().error("Analysis: {},  problem creating Het/Hom ratio", ident());
      return false;

    }

  } else if (unphased_population_){

    if (not syntheticInbreeding()) {

      ExecEnv::log().error("Analysis: {},  problem with synthetic inbreeding analysis", ident());
      return false;

    }

  } else {

    ExecEnv::log().info("Failed to create Filtered Population and Joined Populations Id: {}", ident());

  }

  // Optional check functions.
  //    joinPopulations();
  //    checkPED();

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::MutationAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}


bool kgl::MutationAnalysis::getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters) {

  // Get the output filename
  auto result = named_parameters.find(OUTPUT_FILE_);
  if (result == named_parameters.end()) {
    ExecEnv::log().error("Analytic: {}; Expected Parameter: {} to be defined. {} is deactivated. Available named Parameters:", ident(), OUTPUT_FILE_, ident());
    for (auto const& [parameter_ident, parameter_value] : named_parameters) {

      ExecEnv::log().info("Analysis: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

    }
    return false;
  }
  output_file_name_ = Utility::filePath(result->second, work_directory);

  ExecEnv::log().info("Analysis: {}, initialized with output file: {}", ident(), output_file_name_);

  return true;

}



std::tuple<bool, double> kgl::MutationAnalysis::processFloatField( const std::shared_ptr<const Variant>& variant_ptr
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

      ExecEnv::log().error("MutationAnalysis::processField, Field: {} expected vector size 1, get vector size: {}, vector: {}",
                           field_name, field_vec.size(), vector_str);
      return {false, 0.0};

    }

  } else {

    ExecEnv::log().error("MutationAnalysis::processField, Field: {} not found for Variant: {}",
                         field_name, variant_ptr->output(',',VariantOutputIndex::START_0_BASED, false));

    return {false, 0.0};

  }

}



std::tuple<bool, double> kgl::MutationAnalysis::processStringField(const std::shared_ptr<const Variant>& variant_ptr,
                                                                   const std::string& field_name) const {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, field_name);

  if (field_opt) {

    std::vector<std::string> field_vec = InfoEvidenceAnalysis::varianttoStrings(field_opt.value());

    if (field_vec.size() != 1) {

      ExecEnv::log().error("MutationAnalysis::processStringField, Field: {} expected vector size 1, get vector size: {}",
                           field_name, field_vec.size());
      return {false, 0.0};

    } else {

      try {

        double frequency = std::stod(field_vec.front());
        return {true, frequency };

      }
      catch(std::exception& e) {

        ExecEnv::log().error("MutationAnalysis::processStringField, Field: {} problem converting: {} to double, exception: {}",
                             field_name, field_vec.front(), e.what());
        return {false, 0.0};

      }

    }

  } else {

    ExecEnv::log().error("MutationAnalysis::processField, Field: {} not found for Variant: {}",
                         field_name, variant_ptr->output(',',VariantOutputIndex::START_0_BASED, false));
    return {false, 0.0};
  }

}



std::string kgl::MutationAnalysis::lookupSuperPopulationField(const std::string& super_population) const {

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
void kgl::MutationAnalysis::joinPopulations() const {


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

      ExecEnv::log().warn("MutationAnalysis::checkPED, Genome sample: {} does not have a PED record", genome_id);
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
void kgl::MutationAnalysis::checkPED() const {

  size_t PED_record_count = 0;
  for (auto const& [genome_id, genome_ptr] : diploid_population_->getMap()) {

    auto result = ped_data_->getMap().find(genome_id);

    if (result == ped_data_->getMap().end()) {

      ExecEnv::log().warn("MutationAnalysis::checkPED, Genome sample: {} does not have a PED record", genome_id);

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



