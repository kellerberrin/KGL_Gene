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
#include "kgl_phylogenetic_gene.h"

namespace kgl = kellerberrin::genome;


std::shared_ptr<const kgl::GenomeVariant> getGenomeVariants(std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                                            const std::string& file_name,
                                                            bool vcf_is_gatk,
                                                            const std::string& genome_name,
                                                            kgl::Phred_t read_quality,
                                                            kgl::Phred_t variant_quality,
                                                            long min_count,
                                                            double min_proportion) {

  // Read in the SAM/VCF file variants
  std::shared_ptr<const kgl::GenomeVariant> all_variant_ptr = kgl::VariantFactory().createVariants(genome_db_ptr,
                                                                                                   genome_name,
                                                                                                   file_name,
                                                                                                   vcf_is_gatk,
                                                                                                   read_quality,
                                                                                                   variant_quality,
                                                                                                   min_count,
                                                                                                   min_proportion);

  // Filter on  quality >= 5.
  std::shared_ptr<const kgl::GenomeVariant> filter_ptr = all_variant_ptr->filterVariants(kgl::QualityFilter(5));

  kgl::ExecEnv::log().info("Filtered for quality: {}, Genome: {} has: {} variants", 5, genome_name, filter_ptr->size());

  // Return the genome variants.
  return filter_ptr;

}


kgl::PhylogeneticExecEnv::Application::Application(kgl::Logger& log, const kgl::Phylogenetic& args) {

  log.SetVerbose(args.verbose);

  // Create a genome database object.
  std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr = GenomeDatabase::createGenomeDatabase(args.fastaFile,
                                                                                                  args.gffFile,
                                                                                                  args.gafFile,
                                                                                                  args.aminoTranslationTable);

  // Create a population object.
  std::shared_ptr<kgl::PopulationVariant> pop_variant_ptr(std::make_shared<kgl::PopulationVariant>("Falciparum"));

  // For all organisms
  for (const auto& file : args.fileList) {

    // Generate organism genome variants.
    std::shared_ptr<const kgl::GenomeVariant> variant_ptr = getGenomeVariants(genome_db_ptr,
                                                                              file.file_name,
                                                                              args.vcfAllGATK,
                                                                              file.genome_name,
                                                                              args.readQuality,
                                                                              args.variantQuality,
                                                                              args.minCount,
                                                                              args.minProportion);

    // Store the organism variants in the population object.
    pop_variant_ptr->addGenomeVariant(variant_ptr);



  }

  ApplicationAnalysis::outputSequenceCSV(args.outCSVFile, genome_db_ptr, pop_variant_ptr);

//  GeneAnalysis::mutateAllRegions(args.outCSVFile, 1000,  pop_variant_ptr, genome_db_ptr);

  std::string fasta_file = Utility::filePath(ACTIVE_SEQUENCE, args.workDirectory) + ".fasta";
//  GeneAnalysis::mutateGene(ACTIVE_CONTIG, ACTIVE_GENE, ACTIVE_SEQUENCE, pop_variant_ptr, genome_db_ptr, fasta_file);

//  GeneAnalysis::mutateGenomeRegion("cambodia_fb_SRR2317584", "Pf3D7_02_v3", 534345, 30, pop_variant_ptr, genome_db_ptr);

  // Perform population analysis
  std::string newick_file = Utility::filePath("UPGMA_newick", args.workDirectory) + ".txt";
//  kgl::PhylogeneticAnalysis::UPGMA(newick_file, pop_variant_ptr, genome_db_ptr);

}


