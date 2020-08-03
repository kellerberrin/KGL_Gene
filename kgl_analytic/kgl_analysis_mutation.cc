//
// Created by kellerberrin on 3/7/20.
//

#include "kel_thread_pool.h"
#include "kgl_analysis_mutation.h"
#include "kgl_filter.h"

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

      if (not hetHomRatio(diploid_population)) {

        ExecEnv::log().error("Analysis: {},  problem creating Het/Hom ratio", ident());
        return false;

      }

      filtered_population_ = diploid_population;

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a Diploid Population", ident(), data_object_ptr->Id());
      return false;

    }

  }

  if (data_object_ptr->dataType() == DataTypeEnum::UnphasedPopulation) {

    std::shared_ptr<const UnphasedPopulation> unphased_population = std::dynamic_pointer_cast<const UnphasedPopulation>(data_object_ptr);

    if (unphased_population) {

      filtered_joining_population_ = unphased_population;

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

  if (filtered_population_ and filtered_joining_population_ and ped_data_) {

    ExecEnv::log().info("Filtered Population : {} Count: {} and Joined Population: {} Count: {} both active",
                        filtered_population_->variantCount(), filtered_population_->populationId(),
                        filtered_joining_population_->populationId(), filtered_joining_population_->variantCount());

    joinPopulations();
    checkPED();

  } else {

    ExecEnv::log().info("Failed to create Filtered Population and Joined Populations Id: {}", ident());

  }

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


// Calculate the HetHom Ratio
bool kgl::MutationAnalysis::hetHomRatio(std::shared_ptr<const DiploidPopulation> population) {

  for (auto const& [genome_contig_id, genome_contig_ptr] : genome_GRCh38_->getMap()) {

    ThreadPool thread_pool;
    std::vector<std::future<ret_tuple>> future_vector;
    std::map<GenomeId_t, std::pair<size_t, size_t>> genome_variant_count;
    for (auto const& [genome, genome_ptr] : population->getMap()) {

      std::future<ret_tuple> future = thread_pool.enqueueTask(&MutationAnalysis::processContig, this, genome_contig_id, genome_ptr);
      future_vector.push_back(std::move(future));

    }

    for (auto& future : future_vector) {

      auto [genome, contig_flag, hetero_count, homo_count] = future.get();

      if (contig_flag) {

        genome_variant_count[genome] = std::pair<size_t, size_t>(hetero_count, homo_count);

      }

    }

    if (not genome_variant_count.empty()) {

      // Append the results.
      std::ofstream outfile;
      outfile.open(output_file_name_, std::ofstream::out |  std::ofstream::app);

      outfile << genome_contig_ptr->contigId() << DELIMITER_;
      outfile << genome_variant_count.size() << DELIMITER_;

      for (auto const& [genome_id, het_hom_count] : genome_variant_count) {

        outfile << genome_id << DELIMITER_;

      }

      outfile << '\n';

      outfile << "Hetero" << DELIMITER_;
      outfile << genome_variant_count.size() << DELIMITER_;

      for (auto const& [genome_id, het_hom_count]  : genome_variant_count) {

        outfile << het_hom_count.first << DELIMITER_;

      }

      outfile << '\n';

      outfile << "Homo" << DELIMITER_;
      outfile << genome_variant_count.size() << DELIMITER_;

      for (auto const& [genome_id, het_hom_count] : genome_variant_count) {

        outfile << het_hom_count.second << DELIMITER_;

      }

      outfile << '\n';

      outfile.flush();

    }

  }

  return true;

}


std::tuple<kgl::GenomeId_t, bool, size_t, size_t>
kgl::MutationAnalysis::processContig(ContigId_t contig_id, std::shared_ptr<const DiploidGenome> genome_ptr) {

  auto contig_opt = genome_ptr->getContig(contig_id);

  if (contig_opt) {

    auto contig_ptr = contig_opt.value()->filterVariants(SNPFilter());

    size_t heterozygous_count{0};
    size_t homozygous_count{0};
    for (auto const&[offset, offset_ptr] : contig_ptr->getMap()) {

      OffsetVariantArray variants = offset_ptr->getVariantArray();
      if (variants.size() == 1) {

        ++heterozygous_count;

      } else if (variants.size() == 2) {

        if (variants[0]->homozygous(*variants[1])) {

          ++homozygous_count;

        }

      }

    }

    return { genome_ptr->genomeId(), true, heterozygous_count, homozygous_count};

  }

  return { genome_ptr->genomeId(), false, 0, 0}; // contig not present.

}

// Join Diploid and a single genome population such as Gnomad or Clinvar.
void kgl::MutationAnalysis::joinPopulations() {


  JoinSingleGenome join_pop(filtered_population_, filtered_joining_population_, ped_data_);

  join_pop.joinPopulations();

  ExecEnv::log().info("Joined: {} Genomes in Population: {}", filtered_population_->getMap().size(), filtered_population_->populationId());

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

  class LocalGenomeJoin {

  public :

    explicit LocalGenomeJoin(std::shared_ptr<const GenomeVariant> unphased_genome_ptr) : unphased_genome_ptr_(std::move(unphased_genome_ptr)) {}

    size_t variants_processed_{0};
    size_t joined_variants_found_{0};
    std::shared_ptr<const GenomeVariant> unphased_genome_ptr_;

    // Joins a single genome population (Gnomad, Clinvar) to another (generally phased Diploid) population.
    bool lookupJoinedPop(std::shared_ptr<const Variant> variant_ptr) {

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

  };

  LocalGenomeJoin genome_join(unphased_genome_ptr);

  diploid_genome->processAll(genome_join, &LocalGenomeJoin::lookupJoinedPop);

  return {diploid_genome->genomeId(), genome_join.variants_processed_, genome_join.joined_variants_found_};

}

// Check that all samples (genomes) have a corresponding PED record.
void kgl::MutationAnalysis::checkPED() {

  size_t PED_record_count = 0;
  for (auto const& [genome_id, genome_ptr] : filtered_population_->getMap()) {

    auto result = ped_data_->getMap().find(genome_id);

    if (result == ped_data_->getMap().end()) {

      ExecEnv::log().warn("MutationAnalysis::checkPED, Genome sample: {} does not have a PED record", genome_id);

    } else {

      ++PED_record_count;

    }

  }

  ExecEnv::log().info("Genome samples with PED records: {}", PED_record_count);

}

