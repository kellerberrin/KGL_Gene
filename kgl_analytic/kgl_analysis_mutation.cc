//
// Created by kellerberrin on 3/7/20.
//

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

  };

  if (not getParameters(work_directory, named_parameters)) {

    return false;

  }

// Clear the data file a
  std::ofstream outfile;
  outfile.open(output_file_name_, std::ofstream::out | std::ofstream::trunc);

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::fileReadAnalysis(std::shared_ptr<const PopulationBase> population_base) {

  ExecEnv::log().info("Analysis: {}, begin processing VCF file", ident(), population_base->populationId());

  // Superclass the population
  std::shared_ptr<const DiploidPopulation> diploid_population = std::dynamic_pointer_cast<const DiploidPopulation>(population_base);
  std::shared_ptr<const UnphasedPopulation> unphased_population = std::dynamic_pointer_cast<const UnphasedPopulation>(population_base);

  if (diploid_population) {

    ExecEnv::log().info("Analysis: {}, Generate Hetreozygous/Homozygous ratio statistics, ident()", population_base->populationId());

    if (not hetHomRatio(diploid_population)) {

      ExecEnv::log().error("Analysis: {},  problem creating Het/Hom ratio", ident());
      return false;

    }

//    filtered_population_ = diploid_population->filterVariants(RegionFilter(start_region_, end_region_));
    filtered_population_ = diploid_population;


  }

  if (unphased_population) {

//    filtered_joining_population_ = unphased_population->filterVariants(RegionFilter(start_region_, end_region_));
    filtered_joining_population_ = unphased_population;

  }


  ExecEnv::log().info("Analysis: {}, completed VCF file", ident(), population_base->populationId());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  if (filtered_population_ and filtered_joining_population_) {

    ExecEnv::log().info("Filtered Population : {} Count: {} and Joined Population: {} Count: {} both active",
                        filtered_population_->variantCount(), filtered_population_->populationId(),
                        filtered_joining_population_->populationId(), filtered_joining_population_->variantCount());

    joinPopulations();

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

  for (auto const& [contig, contig_ptr] : genome_GRCh38_->getMap()) {

    std::map<GenomeId_t, std::pair<size_t, size_t>> genome_variant_count;
    for (auto const& [genome, genome_ptr] : population->getMap()) {

      auto contig_opt = genome_ptr->getContig(contig);

      if (contig_opt) {

        auto contig_ptr = contig_opt.value()->filterVariants(SNPFilter());

        std::pair<size_t, size_t> variant_count{0, 0};
        for (auto const& offset_ptr : contig_ptr->getMap()) {

          OffsetVariantArray variants = offset_ptr.second->getVariantArray();
          if (variants.size() == 1) {

            ++variant_count.first;

          } else if (variants.size() == 2) {

            if (variants[0]->homozygous(*variants[1])) {

              ++variant_count.second;

            }

          }

        }

        genome_variant_count[genome] = variant_count;

      }

    }

    if (not genome_variant_count.empty()) {

      // Append the results.
      std::ofstream outfile;
      outfile.open(output_file_name_, std::ofstream::out |  std::ofstream::app);

      outfile << contig_ptr->contigId() << DELIMITER_;
      outfile << genome_variant_count.size() << DELIMITER_;

      for (auto const& genome : genome_variant_count) {

        outfile << genome.first << DELIMITER_;

      }

      outfile << '\n';

      outfile << contig_ptr->contigId() << DELIMITER_;
      outfile << genome_variant_count.size() << DELIMITER_;

      for (auto const& count : genome_variant_count) {

        outfile << count.second.first << DELIMITER_;

      }

      outfile << '\n';

      outfile << contig_ptr->contigId() << DELIMITER_;
      outfile << genome_variant_count.size() << DELIMITER_;

      for (auto const& count : genome_variant_count) {

        outfile << count.second.second << DELIMITER_;

      }

      outfile << '\n';

      outfile.flush();

    }

  }

  return true;

}


// Join Diploid and a signle genome population such as Gnomad or Clinvar.
void kgl::MutationAnalysis::joinPopulations() {


  JoinSingleGenome join_pop(filtered_population_, filtered_joining_population_);

  join_pop.joinPopulations();

  ExecEnv::log().info("Processed Variants: {},  found joined variants: {}",
                      join_pop.variantsProcessed(), join_pop.joinedVariantsFound());

}


bool kgl::JoinSingleGenome::joinPopulations() {

  for (auto const& [genome, genome_ptr] : joined_population_->getMap()) {

    variants_processed_ = 0;
    joined_variants_found_ = 0;

    genome_ptr->processAll(*this, &JoinSingleGenome::lookupJoinedPop);

    ExecEnv::log().info("Genome: {}, Variants: {} Clinvar Variant: {}", genome, variants_processed_, joined_variants_found_);

  }

  return true;

}



// Joins a single genome population (Gnomad, Clinvar) to another (generally phased Diploid) population.
bool kgl::JoinSingleGenome::lookupJoinedPop(std::shared_ptr<const Variant> variant_ptr) {

  ++variants_processed_;

  auto iter = joining_population_->getMap().begin();

  if (iter != joining_population_->getMap().end()) {

    auto contig_opt =  (*iter).second->getContig(variant_ptr->contigId());

    if (contig_opt) {

      auto variant_opt = contig_opt.value()->findVariant(*variant_ptr);

      if (variant_opt) {

        ++joined_variants_found_;

      }

    }

  }

  return true;

}
