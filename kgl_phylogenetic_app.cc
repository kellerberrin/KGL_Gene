//
// Created by kellerberrin on 10/11/17.
//

#include <sstream>
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


std::shared_ptr<const kgl::GenomeVariant> getSNPVariants(kgl::Logger& log,
                                                         std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                                         std::shared_ptr<kgl::PopulationStatistics> population_stats_ptr,
                                                         const std::string& file_name,
                                                         const std::string& genome_name,
                                                         unsigned char read_quality,
                                                         long min_count,
                                                         double min_proportion,
                                                         const std::string& workDirectory) {

  // Read in the SAM file.
  std::shared_ptr<const kgl::ContigCountData> count_data_ptr = kgl::SamCountReader(log).readSAMFile(genome_db_ptr,
                                                                                                    file_name,
                                                                                                    read_quality);

  // Generate all variants.

  std::shared_ptr<const kgl::GenomeVariant> all_variant_ptr = kgl::VariantFactory().createVariants(genome_name,
                                                                                                   count_data_ptr,
                                                                                                   genome_db_ptr,
                                                                                                   min_count,
                                                                                                   min_proportion);

  // Not interested in Synonymous SNPs, so filter them.
  all_variant_ptr = all_variant_ptr->filterVariants(kgl::SynonymousSNPFilter());


  // Create a genome statistics object.
  std::shared_ptr<const kgl::GenomeStatistics> statistics_ptr(std::make_shared<kgl::GenomeStatistics>(genome_db_ptr,
                                                                                                      all_variant_ptr));

  // Write the genome stats to file.
  std::string stats_file_name = kgl::ExecEnv::filePath("genome_stats", workDirectory) + ".csv";
  statistics_ptr->outputCSV(stats_file_name, kgl::VariantOutputIndex::START_1_BASED);

//#define TEST_GENE "PF3D7_1211900"
//#define TEST_GENE "PF3D7_1255200"
//#define TEST_GENE "PF3D7_0101900"
#define TEST_GENE_1 "PF3D7_0711700"   // possible duplicate mutations.
#define TEST_GENE  "PF3D7_0712600"   // possible duplicate mutations.

  // Filter on dup. GENE.
  std::shared_ptr<const kgl::GenomeVariant> filter_ptr = all_variant_ptr->filterVariants(kgl::GeneFilter(TEST_GENE_1));
  filter_ptr = filter_ptr->filterVariants(kgl::CodingFilter());  // remove introns.

  kgl::ExecEnv::log().info("Filtered for: {}, Genome has: {} variants", TEST_GENE_1, filter_ptr->size());
  std::cout << *filter_ptr;

  // Create stats for the GENE.
  std::shared_ptr<const kgl::GenomeStatistics> filtered_stats_ptr(std::make_shared<kgl::GenomeStatistics>(genome_db_ptr,
                                                                                                          filter_ptr));
  // Add to the population stats for later phylogenetic analysis.
  population_stats_ptr->addGenomeStatistics(filtered_stats_ptr);

  // Filter on dup. GENE.
  all_variant_ptr = all_variant_ptr->filterVariants(kgl::GeneFilter(TEST_GENE));
  all_variant_ptr = all_variant_ptr->filterVariants(kgl::CodingFilter());  // remove introns.

  kgl::ExecEnv::log().info("Filtered for: {}, Genome has: {} variants", TEST_GENE, all_variant_ptr->size());
  std::cout << *all_variant_ptr;


  // Return the filtered variants.
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

  for (const auto& file : args.fileList) {

    // Generate mutant filtered simple SNPs.
    std::shared_ptr<const kgl::GenomeVariant> variant_ptr = getSNPVariants(log,
                                                                           genome_db_ptr,
                                                                           population_stats_ptr,
                                                                           file.file_name,
                                                                           file.genome_name,
                                                                           args.readQuality,
                                                                           args.minCount,
                                                                           args.minProportion,
                                                                           args.workDirectory);


    variant_ptr->outputCSV(args.outCSVFile, VariantOutputIndex::START_1_BASED, false);
    std::string fasta_file_name = ExecEnv::filePath(file.genome_name, args.workDirectory) + ".fasta";
    std::string sequence_name = "PF3D7_1211900_" + file.genome_name;
    std::string read_fasta_name = ExecEnv::filePath("PF3D7_1211900_Reference", args.workDirectory) + ".fasta";


//#define MUTANT_PROTEIN 1
//#define MALAWI 1

#ifdef MUTANT_PROTEIN
#ifdef MALAWI

    ApplicationAnalysis::writeMutantProtein(fasta_file_name,
                                            sequence_name,
                                            "Pf3D7_12_v3",
                                            "PF3D7_1211900",
                                            "PF3D7_1211900.1",
                                            genome_db_ptr,
                                            variant_ptr);

    std::string comparison;
    if (ApplicationAnalysis::readMutantProtein(read_fasta_name,
                                           "PF3D7_1211900",
                                           "Pf3D7_12_v3",
                                           "PF3D7_1211900",
                                           "PF3D7_1211900.1",
                                           genome_db_ptr,
                                           variant_ptr,
                                           comparison)) {

      ExecEnv::log().info("PF3D7_1211900 comparison:\n{}", comparison);

    }

#else

    ApplicationAnalysis::writeMutantProtein(fasta_file_name,
                                            sequence_name,
                                            "chr12",
                                            "PF3D7_1211900",
                                            "rna_PF3D7_1211900-1",
                                            genome_db_ptr,
                                            variant_ptr);
    std::string comparison;
    if (ApplicationAnalysis::readMutantProtein(read_fasta_name,
                                           "PF3D7_1211900",
                                           "chr12",
                                           "PF3D7_1211900",
                                           "rna_PF3D7_1211900-1",
                                           genome_db_ptr,
                                           variant_ptr,
                                           comparison)) {

      ExecEnv::log().info("PF3D7_1211900 comparison:\n{}", comparison);

    }

#endif
#endif

    // Generate a UPGMA newick file from the population stats object
    std::string newick_file = ExecEnv::filePath("UPGMA_newick", args.workDirectory) + ".txt";
    kgl::PhylogeneticAnalysis::UPGMA(newick_file, population_stats_ptr);

  } // For all sam files loop.

  // Perform population analysis.


}


