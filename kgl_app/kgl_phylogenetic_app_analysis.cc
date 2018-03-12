//
// Created by kellerberrin on 12/02/18.
//


#include <kgl_sequence_distance.h>
#include "kgl_phylogenetic_app.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_phylogenetic_gene.h"
#include "kgl_upgma.h"
#include "kgl_rna_search.h"
#include "kgl_ploidy_analysis.h"


namespace kgl = kellerberrin::genome;



void kgl::PhylogeneticApp::performAnalysis(const kgl::Phylogenetic& args,
                                           std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                           std::shared_ptr<const kgl::PopulationVariant> pop_variant_ptr) {


  if (args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("No analytic specified - performing all analytics (time consuming)");

  }

  if (args.analysisType == kgl::Phylogenetic::ANALYZE_SEQUENCES or args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("Analyzing coding sequences");
    std::shared_ptr<const GlobalAminoSequenceDistance> amino_distance_metric(std::make_shared<const LevenshteinGlobal>());
    std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric(std::make_shared<const LevenshteinGlobal>());
    std::string coding_file = kgl::Utility::filePath("CodingAnalysis", args.workDirectory) + ".csv";
    kgl::ApplicationAnalysis::outputSequenceCSV(coding_file, dna_distance_metric, amino_distance_metric, genome_db_ptr, pop_variant_ptr);

  }

  if (args.analysisType == kgl::Phylogenetic::ANALYZE_INTERVAL or args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("Analyzing genome intervals");
    std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric(std::make_shared<const LevenshteinGlobal>());
    std::string interval_file = kgl::Utility::filePath("IntervalAnalysis", args.workDirectory) + ".csv";
    kgl::GeneAnalysis::mutateAllRegions(interval_file, 1000,  dna_distance_metric, pop_variant_ptr, genome_db_ptr);

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

  std::string DNA_mutation_file = kgl::Utility::filePath("DNAMutations", args.workDirectory) + ".csv";
  ApplicationAnalysis::outputDNAMutationCSV(DNA_mutation_file, PFATP4_CONTIG, PFATP4_GENE, PFATP4_SEQUENCE, genome_db_ptr, pop_variant_ptr);
  std::string amino_mutation_file = kgl::Utility::filePath("AminoMutations", args.workDirectory) + ".csv";
  ApplicationAnalysis::outputAminoMutationCSV(amino_mutation_file, PFATP4_CONTIG, PFATP4_GENE, PFATP4_SEQUENCE, genome_db_ptr, pop_variant_ptr);

  if (args.analysisType == kgl::Phylogenetic::ANALYZE_UPGMA or args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("Performing a UPGMA analytic");
    std::string newick_file = "UPGMA_newick.txt";
    std::shared_ptr<const AminoSequenceDistance> distance_metric_ptr(std::make_shared<const Blosum80Global>());
    kgl::UPGMAGenePhyloTree<kgl::UPGMAATP4Distance>(args.workDirectory,
                                                         newick_file,
                                                         distance_metric_ptr,
                                                         pop_variant_ptr,
                                                         genome_db_ptr,
                                                         kgl::UPGMAProteinDistance::SYMBOLIC_ATP4_FAMILY);
  }


  if (args.analysisType == kgl::Phylogenetic::ANALYZE_RNA or args.analysisType == kgl::Phylogenetic::WILDCARD) {

    kgl::ExecEnv::log().info("Performing an RNA analytic");

    RNAAnalysis rna_analysis;

    rna_analysis.getRNARegions("Pf3D7_04_v3",
                               956848,
                               135,
                               StrandSense::FORWARD,
                               "Pf3D7_04_v3",
                               (958066 - 1000),
                               9000,
                               StrandSense::FORWARD,
                               genome_db_ptr);

    std::shared_ptr<const LocalDNASequenceCompare>  compare_metric_ptr(std::make_shared<const DNALocalAffineGap>());
    rna_analysis.compareRNARegion(0, 24, 1, compare_metric_ptr);

    rna_analysis.showResults(2);

  }

}



