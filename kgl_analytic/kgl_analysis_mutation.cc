//
// Created by kellerberrin on 3/7/20.
//

#include "kgl_analysis_mutation.h"
#include "kgl_variant_db_phased.h"
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

  // Superclass the population
  std::shared_ptr<const DiploidPopulation> population = std::dynamic_pointer_cast<const DiploidPopulation>(population_base);

  if (not population) {

    ExecEnv::log().error("Analysis: {},  expected a Diploid Population", ident());
    return false;

  }

  for (auto const& [contig, contig_ptr] : genome_GRCh38_->getMap()) {

    std::map<GenomeId_t, std::pair<size_t, size_t>> genome_variant_count;
    for (auto const& [genome, genome_ptr] : population->getMap()) {

      auto contig_opt = genome_ptr->getContig(contig);

      if (contig_opt) {

        auto contig_ptr = contig_opt.value()->filterVariants(SNPFilter());

        std::pair<size_t, size_t> variant_count{0, 0};
        for (auto const& [offset, offset_ptr] : contig_ptr->getMap()) {

          OffsetVariantArray variants = offset_ptr->getVariantArray();
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

      for (auto const& [genome, count] : genome_variant_count) {

        outfile << genome << DELIMITER_;

      }

      outfile << '\n';

      outfile << contig_ptr->contigId() << DELIMITER_;
      outfile << genome_variant_count.size() << DELIMITER_;

      for (auto const& [genome, count] : genome_variant_count) {

        outfile << count.first << DELIMITER_;

      }

      outfile << '\n';

      outfile << contig_ptr->contigId() << DELIMITER_;
      outfile << genome_variant_count.size() << DELIMITER_;

      for (auto const& [genome, count] : genome_variant_count) {

        outfile << count.second << DELIMITER_;

      }

      outfile << '\n';

      outfile.flush();

    }

  }

  ExecEnv::log().info("Analysis: {}, completed VCF file", ident(), population->populationId());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

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
