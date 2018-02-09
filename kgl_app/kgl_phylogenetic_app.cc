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


void performAnalysis(const kgl::Phylogenetic& args,
                     std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                     std::shared_ptr<const kgl::PopulationVariant> pop_variant_ptr) {


  if (args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("No analytic specified - performing all analytics (time consuming)");

  }

  if (args.analysisType == kgl::Phylogenetic::ANALYZE_SEQUENCES or args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("Analyzing coding sequences");
    kgl::ApplicationAnalysis::outputSequenceCSV(args.outCSVFile, genome_db_ptr, pop_variant_ptr);

  }

  if (args.analysisType == kgl::Phylogenetic::ANALYZE_INTERVAL or args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("Analyzing genome intervals");
    kgl::GeneAnalysis::mutateAllRegions(args.outCSVFile, 1000,  pop_variant_ptr, genome_db_ptr);

  }

  if (args.analysisType == kgl::Phylogenetic::ANALYZE_GENE or args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("Analyzing a specific sequence");
    std::string fasta_file = kgl::Utility::filePath(ACTIVE_SEQUENCE, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGene(ACTIVE_CONTIG, ACTIVE_GENE, ACTIVE_SEQUENCE, pop_variant_ptr, genome_db_ptr, fasta_file);

  }


  if (args.analysisType == kgl::Phylogenetic::ANALYZE_REGION or args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("Analyzing genome regions");

    std::string region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_944950";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 944950, 135, pop_variant_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_941875";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 941875, 4293, pop_variant_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_569342";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 569342, 7467, pop_variant_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_584668";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 584668, 7280, pop_variant_ptr, genome_db_ptr, region_fasta_file);

  }

  if (args.analysisType == kgl::Phylogenetic::ANALYZE_UPGMA or args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("Performing a UPGMA analytic");
    // Perform population analysis
    std::string newick_file = kgl::Utility::filePath("UPGMA_newick", args.workDirectory) + ".txt";
    kgl::PhylogeneticAnalysis::UPGMA(newick_file, pop_variant_ptr, genome_db_ptr);

  }

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

  performAnalysis(args, genome_db_ptr, pop_variant_ptr);



}


