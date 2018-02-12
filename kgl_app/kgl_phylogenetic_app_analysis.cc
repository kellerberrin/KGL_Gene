//
// Created by kellerberrin on 12/02/18.
//


#include "kgl_phylogenetic_app.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_phylogenetic_gene.h"
#include "kgl_upgma.h"
#include "kgl_utility.h"


namespace kgl = kellerberrin::genome;



void kgl::PhylogeneticApp::performAnalysis(const kgl::Phylogenetic& args,
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
    region_fasta_file += "_561666";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 561666, 11753, pop_variant_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_462000";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_01_v3", 462000, 2000, pop_variant_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_0";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 0, 1200490, pop_variant_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
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
    std::string newick_file = kgl::Utility::filePath("UPGMA_newick", args.workDirectory) + ".txt";
    kgl::UPGMATree<kgl::UPGMAProteinDistance, std::string>(newick_file,
                                                          pop_variant_ptr,
                                                          genome_db_ptr,
                                                          kgl::UPGMAProteinDistance::SYMBOLIC_RIFIN_FAMILY);
  }

}


