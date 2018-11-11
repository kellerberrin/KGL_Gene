//
// Created by kellerberrin on 12/02/18.
//


#include <kgl_read_phasing.h>
#include "kgl_sequence_distance.h"
#include "kgl_genome_aux_csv.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_phylogenetic_gene.h"
#include "kgl_upgma.h"
#include "kgl_rna_search.h"
#include "kgl_phylogenetic_app_analysis.h"
#include "kgl_variant_classify.h"
#include "kgl_variant_phase.h"
#include "kgl_distribution_analysis.h"
#include "kgl_finestructure_analysis.h"
#include "kgl_upgma_unphased.h"
#include "kgl_epigenetic_motif.h"

#include "kgd_deconvolv_app.h"


namespace kgl = kellerberrin::genome;
namespace kgd = kellerberrin::deconvolv;


bool kgl::PhylogeneticAnalysis::checkAnalysisType(const std::string& analysis_type) {


  if (analysis_type != ANALYZE_INTERVAL
      and analysis_type != ANALYZE_SEQUENCES
      and analysis_type != ANALYZE_GENE
      and analysis_type != ANALYZE_REGION
      and analysis_type != ANALYZE_UPGMA
      and analysis_type != ANALYZE_RNA
      and analysis_type != ANALYZE_SNP
      and analysis_type != ANALYZE_MIX
      and analysis_type != ANALYZE_MOTIF) {

    ExecEnv::log().error("Invalid Analysis Type: {}.  Must be one of: {}, {}, {}, {}, {}, {}, {}, {}, {}.",
                         analysis_type,
                         ANALYZE_INTERVAL,
                         ANALYZE_SEQUENCES,
                         ANALYZE_GENE,
                         ANALYZE_REGION,
                         ANALYZE_UPGMA,
                         ANALYZE_RNA,
                         ANALYZE_SNP,
                         ANALYZE_MIX,
                         ANALYZE_MOTIF);

    return false;

  }


  return true;

}




void kgl::PhylogeneticAnalysis::performAnalysis(const kgl::Phylogenetic& args,
                                                const kgl::PropertyTree& runtime_options,
                                                std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                                std::shared_ptr<const UnphasedPopulation> unphased_population_ptr,
                                                std::shared_ptr<const PhasedPopulation> population_ptr) {


  // Split into country populations.
  std::vector<CountryPair> country_pairs = GenomeAuxData::getCountries(args.auxCSVFile, population_ptr);



  if (args.analysisType == ANALYZE_MOTIF) {

    std::string motif_file = kgl::Utility::filePath("Motif", args.workDirectory) + ".txt";
    PromoterMotif::displayTFFMotif(genome_db_ptr, motif_file);

  }

  if (args.analysisType == ANALYZE_SEQUENCES) {

    kgl::ExecEnv::log().info("Analyzing coding sequences");

    std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric(std::make_shared<const LevenshteinGlobal>());

    std::string coding_file = kgl::Utility::filePath("All_DNA_CodingAnalysis", args.workDirectory) + ".csv";
    kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::DNA, dna_distance_metric,  genome_db_ptr, population_ptr);

    coding_file = kgl::Utility::filePath("All_Variant_CodingAnalysis", args.workDirectory) + ".csv";
    kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::VARIANT, dna_distance_metric,  genome_db_ptr, population_ptr);

    coding_file = kgl::Utility::filePath("All_SNP_CodingAnalysis", args.workDirectory) + ".csv";
    kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SNP, dna_distance_metric,  genome_db_ptr, population_ptr);

    coding_file = kgl::Utility::filePath("All_SIZE_CodingAnalysis", args.workDirectory) + ".csv";
    kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SIZE, dna_distance_metric,  genome_db_ptr, population_ptr);
    
    std::shared_ptr<const GlobalAminoSequenceDistance> amino_distance_metric(std::make_shared<const LevenshteinGlobal>());
    coding_file = kgl::Utility::filePath("All_Amino_CodingAnalysis", args.workDirectory) + ".csv";
    kgl::ApplicationAnalysis::outputAminoSequenceCSV(coding_file, amino_distance_metric, genome_db_ptr, population_ptr);

    for (auto country : country_pairs) {

      coding_file = kgl::Utility::filePath(country.first + "Variant_CodingAnalysis", args.workDirectory) + ".csv";
      kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::VARIANT, dna_distance_metric,  genome_db_ptr, country.second);

      coding_file = kgl::Utility::filePath(country.first + "SNP_CodingAnalysis", args.workDirectory) + ".csv";
      kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SNP, dna_distance_metric,  genome_db_ptr, country.second);

      coding_file = kgl::Utility::filePath(country.first + "SIZE_CodingAnalysis", args.workDirectory) + ".csv";
      kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SIZE, dna_distance_metric,  genome_db_ptr, country.second);

      coding_file = kgl::Utility::filePath(country.first + "DNA_CodingAnalysis", args.workDirectory) + ".csv";
      kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::DNA, dna_distance_metric,  genome_db_ptr, country.second);

      coding_file = kgl::Utility::filePath(country.first + "Amino_CodingAnalysis", args.workDirectory) + ".csv";
      kgl::ApplicationAnalysis::outputAminoSequenceCSV(coding_file, amino_distance_metric, genome_db_ptr, country.second);

    }

  }

  if (args.analysisType == ANALYZE_INTERVAL) {

#define INTERVAL_SIZE 100

    kgl::ExecEnv::log().info("Analyzing genome intervals");

    AggregateVariantDistribution variant_distribution;
    variant_distribution.variantDistribution(population_ptr);
    std::string file_name = "IntervalAnalysis_all";
    std::string interval_file = kgl::Utility::filePath(file_name, args.workDirectory) + ".csv";
    variant_distribution.writeDistribution(genome_db_ptr, INTERVAL_SIZE, interval_file, ',');

    for (auto country : country_pairs) {

      AggregateVariantDistribution variant_distribution;
      variant_distribution.variantDistribution(country.second);
      std::string file_name = "IntervalAnalysis_" + country.first;
      std::string interval_file = kgl::Utility::filePath(file_name, args.workDirectory) + ".csv";
      variant_distribution.writeDistribution(genome_db_ptr, INTERVAL_SIZE, interval_file, ',');

    }

  }

  if (args.analysisType == ANALYZE_GENE) {

    kgl::ExecEnv::log().info("Analyzing a specific sequence");
    std::string fasta_file = kgl::Utility::filePath(ACTIVE_SEQUENCE, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGene(ACTIVE_CONTIG, ACTIVE_GENE, ACTIVE_SEQUENCE, population_ptr, genome_db_ptr, fasta_file);

  }


  if (args.analysisType == ANALYZE_REGION) {

    kgl::ExecEnv::log().info("Analyzing genome regions");

    std::string region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_561666";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 561666, 11753, population_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_462000";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_01_v3", 462000, 2000, population_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_0";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 0, 1200490, population_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_944950";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 944950, 135, population_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_941875";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 941875, 4293, population_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_569342";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 569342, 7467, population_ptr, genome_db_ptr, region_fasta_file);

    region_fasta_file = "malawi_fb_SRR609075";
    region_fasta_file += "_584668";
    region_fasta_file = kgl::Utility::filePath(region_fasta_file, args.workDirectory) + ".fasta";
    kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 584668, 7280, population_ptr, genome_db_ptr, region_fasta_file);

  }

  if (args.analysisType == ANALYZE_MIX) {

    // Index variants by ref/alt count.
    VariantClassifier classifier(unphased_population_ptr);

    for (auto genome : classifier.getGenomes()) {

      const size_t minimum_offset_base_count = 25;
      const size_t maximum_offset_base_count = 150;
      kgd::MixtureDataObj mixture_data = classifier.convertToMixture(genome, minimum_offset_base_count, maximum_offset_base_count);
      kgd::ExecEnv::executeLib(mixture_data);

    }

  }

  if (args.analysisType == ANALYZE_SNP) {

    GenomeAuxData aux_data;
    aux_data.readParseAuxData(args.auxCSVFile);

    std::string DNA_mutation_file = kgl::Utility::filePath("DNAMutations.csv", args.workDirectory);
    ApplicationAnalysis::outputDNAMutationCSV(DNA_mutation_file,
                                              PFATP4_CONTIG,
                                              PFATP4_GENE,
                                              PFATP4_SEQUENCE,
                                              genome_db_ptr,
                                              population_ptr,
                                              aux_data);
    std::string amino_mutation_file = kgl::Utility::filePath("AminoMutations.csv", args.workDirectory);
    ApplicationAnalysis::outputAminoMutationCSV(amino_mutation_file, PFATP4_CONTIG, PFATP4_GENE, PFATP4_SEQUENCE,
                                                genome_db_ptr, population_ptr);

//    phased_statistics_ptr->outputPopulation();

  }

  if (args.analysisType == ANALYZE_UPGMA) {

    kgl::ExecEnv::log().info("Generating FineStructure Files.");

#define BASES_PER_CENTIMORGAN 15000.0

    for (auto country : country_pairs) {

      std::string file_name = "FS/" + country.first;
      std::string fine_structure_file = kgl::Utility::filePath(file_name, args.workDirectory);
      FineStructureAnalysis::generateFiles(fine_structure_file, country.second, BASES_PER_CENTIMORGAN);

    }


    kgl::ExecEnv::log().info("Performing a UPGMA analytic");
//    std::shared_ptr<const AminoSequenceDistance> distance_metric_ptr(std::make_shared<const Blosum80Global>());
//    kgl::UPGMAGenePhyloTree<kgl::UPGMAATP4Distance>(args.workDirectory,
//                                                    newick_file,
//                                                    distance_metric_ptr,
//                                                    population_ptr,
//                                                    genome_db_ptr,
//                                                    kgl::UPGMAProteinDistance::SYMBOLIC_ATP4_FAMILY);

//    std::string newick_file = kgl::Utility::filePath("newick", args.workDirectory) + ".txt";
//    std::shared_ptr<const UnphasedPopulation> filtered_ptr = unphased_population_ptr->filterVariants(DPCountFilter(25));
//    kgl::UPGMAUnphasedTree<kgl::UPGMAUnphasedDistance>(newick_file, filtered_ptr);

  }


  if (args.analysisType == ANALYZE_RNA) {

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


