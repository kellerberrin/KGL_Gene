//
// Created by kellerberrin on 10/11/17.
//

#include <sstream>
#include "kgl_utility.h"
#include "kgl_genome_types.h"
#include "kgl_phylogenetic_env.h"
#include "kgl_genome_db.h"
#include "kgl_gff_fasta.h"
#include "kgl_sam_process.h"
#include "kgl_variant_factory.h"
#include "kgl_variant_factory_compound.h"
#include "kgl_filter.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_statistics.h"

namespace kgl = kellerberrin::genome;


std::shared_ptr<const kgl::GenomeVariant> getGenomeVariants(std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                                            std::shared_ptr<kgl::PopulationStatistics> population_stats_ptr,
                                                            const std::string& file_name,
                                                            const std::string& genome_name,
                                                            kgl::Phred_t read_quality,
                                                            kgl::Phred_t variant_quality,
                                                            long min_count,
                                                            double min_proportion,
                                                            const std::string& workDirectory) {

  // Read in the SAM file variants
  std::shared_ptr<const kgl::GenomeVariant> all_variant_ptr = kgl::VariantFactory().createVariants(genome_db_ptr,
                                                                                                   genome_name,
                                                                                                   file_name,
                                                                                                   read_quality,
                                                                                                   variant_quality,
                                                                                                   min_count,
                                                                                                   min_proportion);

  // Create a genome statistics object.
  std::shared_ptr<const kgl::GenomeStatistics> statistics_ptr(std::make_shared<kgl::GenomeStatistics>(genome_db_ptr,
                                                                                                      all_variant_ptr));
  // Add to the population statistics
  population_stats_ptr->addGenomeStatistics(statistics_ptr);

  // Write the genome stats to file.
  std::string stats_file_name = kgl::Utility::filePath("genome_stats", workDirectory) + ".csv";
  statistics_ptr->outputFeatureCSV(stats_file_name, kgl::VariantOutputIndex::START_1_BASED);

  // Filter on sequence
//  std::shared_ptr<const kgl::GenomeVariant> check_ptr = all_variant_ptr->filterVariants(kgl::AndFilter(kgl::ContigFilter("Pf3D7_01_v3"), kgl::RegionFilter(63100,63280)));
//  kgl::ExecEnv::log().info("Region Filter\n:");
//  std::cout << *check_ptr;

// pfATP4 drug target ATP4 sodium pump.
#define PFATP4_MINORITY_CONTIG "chr12"
#define PFATP4_MINORITY_GENE "PF3D7_1211900"
#define PFATP4_MINORITY_SEQUENCE "rna_PF3D7_1211900-1"

#define PFATP4_MALAWI_CONTIG "PF3D7_12_v3"
#define PFATP4_MALAWI_GENE "PF3D7_1211900"
#define PFATP4_MALAWI_SEQUENCE "PF3D7_1211900.1"

// Erythrocyte membrane (Blood Cell Surface) protein (busy Malawi)
#define _BCS_CONTIG "PF3D7_12_v3"
#define _BCS_GENE "PF3D7_1255200"
#define  BCS_SEQUENCE "PF3D7_1255200.1"

// Mutant rich Rifin (Malawi)
#define RIFIN_3_CONTIG "Pf3D7_01_v3"
#define RIFIN_3_GENE "PF3D7_0101900"
#define RIFIN_3_SEQUENCE "PF3D7_0101900.1"

// Mutant rich Rifin (Malawi)
#define RIFIN_4_CONTIG "Pf3D7_08_v3"
#define RIFIN_4_GENE "PF3D7_0808900"
#define RIFIN_4_SEQUENCE "PF3D7_0808900.1"

// Rifin very similar mutations (Malawi).
#define RIFIN_1_CONTIG "Pf3D7_07_v3"
#define RIFIN_1_GENE "PF3D7_0711700"
#define RIFIN_1_SEQUENCE "PF3D7_0711700.1"

#define RIFIN_2_CONTIG "Pf3D7_07_v3"
#define RIFIN_2_GENE  "PF3D7_0712600"
#define RIFIN_2_SEQUENCE "PF3D7_0712600.1"

// S-Antigen very busy 5-prime and 3-prime regions (Malawi).
#define S_ANTIGEN_CONTIG "Pf3D7_10_v3"
#define S_ANTIGEN_GENE "PF3D7_1035200"
#define S_ANTIGEN_SEQUENCE "PF3D7_1035200.1"

#define ACTIVE_CONTIG S_ANTIGEN_CONTIG
#define ACTIVE_GENE S_ANTIGEN_GENE
#define ACTIVE_SEQUENCE S_ANTIGEN_SEQUENCE

  // Filter on sequence
  std::shared_ptr<const kgl::GenomeVariant> filter_ptr = all_variant_ptr->filterVariants(kgl::SequenceFilter(ACTIVE_SEQUENCE));
  // Filter on contig

  kgl::ExecEnv::log().info("Filtered for: {}, Genome: {} has: {} variants", ACTIVE_SEQUENCE, genome_name, filter_ptr->size());
  std::cout << *filter_ptr;

  // Return the genome variants.
  return all_variant_ptr;

}


kgl::PhylogeneticExecEnv::Application::Application(kgl::Logger& log, const kgl::Phylogenetic& args) {

  log.SetVerbose(args.verbose);

  // Create a genome database object.
  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = kgl::ParseGffFasta().readFastaGffFile(args.fastaFile,
                                                                                             args.gffFile);
  // Set the amino translation table
  genome_db_ptr->setTranslationTable(args.aminoTranslationTable);

  // Wire-up the genome database.
  genome_db_ptr->createVerifyGenomeDatabase();

  // Create a population statistics object.
  std::shared_ptr<kgl::PopulationStatistics> population_stats_ptr(std::make_shared<kgl::PopulationStatistics>("Malawi-PRJNA173723"));

  // Create a population variant object.
  std::shared_ptr<kgl::PopulationVariant> population_variant_ptr(std::make_shared<kgl::PopulationVariant>("Malawi-PRJNA173723"));

  // For all organisms
  for (const auto& file : args.fileList) {

    // Generate all genome variants.
    std::shared_ptr<const kgl::GenomeVariant> variant_ptr = getGenomeVariants(genome_db_ptr,
                                                                              population_stats_ptr,
                                                                              file.file_name,
                                                                              file.genome_name,
                                                                              args.readQuality,
                                                                              args.variantQuality,
                                                                              args.minCount,
                                                                              args.minProportion,
                                                                              args.workDirectory);

    // Store the genome variant pointer
    population_variant_ptr->addGenomeVariant(variant_ptr);

    // Write genome variants to file.
    variant_ptr->outputCSV(args.outCSVFile, VariantOutputIndex::START_1_BASED, false);
    std::string sequence_name = file.genome_name;
    sequence_name += "_";
    sequence_name += ACTIVE_GENE;
    std::string fasta_file_name = Utility::filePath(sequence_name, args.workDirectory) + ".fasta";

    // Write the vector of mutant proteins to a fasta file.
    ApplicationAnalysis::writeMutantProteins(fasta_file_name,
                                             sequence_name,
                                             ACTIVE_CONTIG,
                                             ACTIVE_GENE,
                                             ACTIVE_SEQUENCE,
                                             genome_db_ptr,
                                             variant_ptr);

    std::vector<std::string> comparison_vector;

    // Generate a vector of 5 Prime UTR mutation maps for visual inspection
    if (ApplicationAnalysis::compare5Prime(ACTIVE_CONTIG,
                                           ACTIVE_GENE,
                                           ACTIVE_SEQUENCE,
                                           1000,
                                           genome_db_ptr,
                                           variant_ptr,
                                           comparison_vector)) {

      for (const auto& comparison : comparison_vector) {

        ExecEnv::log().info("Genome: {} 5 Prime UTR Sequence:{} Comparison:\n{}\n", file.genome_name, ACTIVE_SEQUENCE, comparison);

      }

    }

    // Generate a vector of protein mutation maps for visual inspection
    if (ApplicationAnalysis::compareMutantProteins(ACTIVE_CONTIG,
                                                   ACTIVE_GENE,
                                                   ACTIVE_SEQUENCE,
                                                   genome_db_ptr,
                                                   variant_ptr,
                                                   comparison_vector)) {

      for (const auto& comparison : comparison_vector) {

        ExecEnv::log().info("Genome: {} Protein Sequence:{} Comparison:\n{}\n", file.genome_name, ACTIVE_SEQUENCE, comparison);

      }

    }

    // Generate a vector of 3 Prime UTR mutation maps for visual inspection
    if (ApplicationAnalysis::compare3Prime(ACTIVE_CONTIG,
                                           ACTIVE_GENE,
                                           ACTIVE_SEQUENCE,
                                           1000,
                                           genome_db_ptr,
                                           variant_ptr,
                                           comparison_vector)) {

      for (const auto& comparison : comparison_vector) {

        ExecEnv::log().info("Genome: {} 3 Prime UTR Sequence:{} Comparison:\n{}\n", file.genome_name, ACTIVE_SEQUENCE, comparison);

      }

    }

    // Generate a vector of 3 Prime UTR mutation maps for visual inspection
    if (ApplicationAnalysis::compareMutantRegions(ACTIVE_CONTIG, 1394830, 20, StrandSense::FORWARD, genome_db_ptr, variant_ptr, comparison_vector)) {

      for (const auto& comparison : comparison_vector) {

        ExecEnv::log().info("Genome: {} Test Sequence:{} Comparison:\n{}\n", file.genome_name, ACTIVE_SEQUENCE, comparison);

      }

    }

    if (ApplicationAnalysis::compareMutantRegions(ACTIVE_CONTIG, 1396590, 20, StrandSense::FORWARD, genome_db_ptr, variant_ptr, comparison_vector)) {

      for (const auto& comparison : comparison_vector) {

        ExecEnv::log().info("Genome: {} Test Sequence:{} Comparison:\n{}\n", file.genome_name, ACTIVE_SEQUENCE, comparison);

      }

    }


  } // For all sam files loop.



  // Perform population analysis
  // (disabled to save CPU time)
  std::string newick_file = Utility::filePath("file.genome_name+ UPGMA_newick", args.workDirectory) + ".txt";
//  kgl::PhylogeneticAnalysis::UPGMA(newick_file, population_stats_ptr);

}


