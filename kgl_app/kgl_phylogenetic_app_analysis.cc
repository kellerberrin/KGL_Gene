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



void kgl::PhylogeneticAnalysis::performAnalysis(const std::string& analysis_type) {

  auto result = analysis_map_.find(analysis_type);

  if (result == analysis_map_.end()) {

    ExecEnv::log().error("Invalid Analysis Type: '{}'. Must be one of the following (case insensitive):", analysis_type);

    for (auto analysis: analysis_map_) {

      ExecEnv::log().info("Valid Analysis Type: '{}'", analysis.first);

    }

    return;

  }

  kgl::ExecEnv::log().info("**********Performing Analysis: {} **************", result->first);

  // Call the analysis.
  (this->*(result->second))();

  kgl::ExecEnv::log().info("**********Completed Analysis: {} **************", result->first);

}



void kgl::PhylogeneticAnalysis::performMotif() {


  std::string motif_file = kgl::Utility::filePath("Motif", runtime_options_.workDirectory()) + ".csv";
  PromoterMotif::displayTFFMotif(get3D7Genome(), motif_file, ',');

}


void kgl::PhylogeneticAnalysis::performSequence() {


  std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric(std::make_shared<const LevenshteinGlobal>());

  std::string coding_file = kgl::Utility::filePath("All_DNA_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::DNA, dna_distance_metric,
                                                 get3D7Genome(), population_ptr_);

  coding_file = kgl::Utility::filePath("All_Variant_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::VARIANT, dna_distance_metric,
                                                 get3D7Genome(), population_ptr_);

  coding_file = kgl::Utility::filePath("All_SNP_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SNP, dna_distance_metric,
                                                 get3D7Genome(), population_ptr_);

  coding_file = kgl::Utility::filePath("All_SIZE_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SIZE, dna_distance_metric,
                                                 get3D7Genome(), population_ptr_);

  std::shared_ptr<const GlobalAminoSequenceDistance> amino_distance_metric(
      std::make_shared<const LevenshteinGlobal>());
  coding_file = kgl::Utility::filePath("All_Amino_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  kgl::ApplicationAnalysis::outputAminoSequenceCSV(coding_file, amino_distance_metric, get3D7Genome(),
                                                   population_ptr_);

  // Split into country populations.
  std::string aux_file_path;
  if (not runtime_options_.getPropertiesAuxFile(aux_file_path)) {

    ExecEnv::log().critical("Aux Genome information file not found");

  }

  std::vector<CountryPair> country_pairs = GenomeAuxData::getCountries(aux_file_path, population_ptr_);

  for (auto country : country_pairs) {

    coding_file = kgl::Utility::filePath(country.first + "Variant_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::VARIANT, dna_distance_metric,
                                                   get3D7Genome(), country.second);

    coding_file = kgl::Utility::filePath(country.first + "SNP_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SNP, dna_distance_metric,
                                                   get3D7Genome(), country.second);

    coding_file = kgl::Utility::filePath(country.first + "SIZE_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SIZE, dna_distance_metric,
                                                   get3D7Genome(), country.second);

    coding_file = kgl::Utility::filePath(country.first + "DNA_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    kgl::ApplicationAnalysis::outputDNASequenceCSV(coding_file, SequenceAnalysisType::DNA, dna_distance_metric,
                                                   get3D7Genome(), country.second);

    coding_file = kgl::Utility::filePath(country.first + "Amino_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    kgl::ApplicationAnalysis::outputAminoSequenceCSV(coding_file, amino_distance_metric, get3D7Genome(),
                                                     country.second);

  }

}

void kgl::PhylogeneticAnalysis::performInterval() {


#define INTERVAL_SIZE 100
#define ANALYSIS_CONTIG "Pf3D7_04_v3"
#define ANALYSIS_START 540000
#define ANALYSIS_END 604000


  AggregateVariantDistribution variant_distribution;
  variant_distribution.variantDistribution(population_ptr_);
  std::string file_name = "IntervalAnalysis_all";
  std::string interval_file = kgl::Utility::filePath(file_name, runtime_options_.workDirectory()) + ".csv";
  variant_distribution.writeDistribution(get3D7Genome(), INTERVAL_SIZE, ANALYSIS_CONTIG, ANALYSIS_START, ANALYSIS_END, true, interval_file, ',');

  // Split into country populations.
  std::string aux_file_path;
  if (not runtime_options_.getPropertiesAuxFile(aux_file_path)) {

    ExecEnv::log().critical("Aux Genome information file not found");

  }

  std::vector<CountryPair> country_pairs = GenomeAuxData::getCountries(aux_file_path, population_ptr_);

  for (auto country : country_pairs) {

    AggregateVariantDistribution variant_distribution;
    variant_distribution.variantDistribution(country.second);
    std::string file_name = "IntervalAnalysis_" + country.first;
    std::string interval_file = kgl::Utility::filePath(file_name, runtime_options_.workDirectory()) + ".csv";
    variant_distribution.writeDistribution(get3D7Genome(), INTERVAL_SIZE, ANALYSIS_CONTIG, 0, 0, true, interval_file, ',');

  }


}


void kgl::PhylogeneticAnalysis::performGene() {

    std::string fasta_file = kgl::Utility::filePath(ACTIVE_SEQUENCE, runtime_options_.workDirectory()) + ".fasta";
    kgl::GeneAnalysis::mutateGene(ACTIVE_CONTIG, ACTIVE_GENE, ACTIVE_SEQUENCE, population_ptr_, get3D7Genome(),
                                  fasta_file);

}


void kgl::PhylogeneticAnalysis::performRegion() {

  std::string region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_561666";
  region_fasta_file = kgl::Utility::filePath(region_fasta_file, runtime_options_.workDirectory()) + ".fasta";
  kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 561666, 11753, population_ptr_,
                                        get3D7Genome(), region_fasta_file);

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_462000";
  region_fasta_file = kgl::Utility::filePath(region_fasta_file, runtime_options_.workDirectory()) + ".fasta";
  kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_01_v3", 462000, 2000, population_ptr_,
                                        get3D7Genome(), region_fasta_file);

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_0";
  region_fasta_file = kgl::Utility::filePath(region_fasta_file, runtime_options_.workDirectory()) + ".fasta";
  kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 0, 1200490, population_ptr_,
                                        get3D7Genome(), region_fasta_file);

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_944950";
  region_fasta_file = kgl::Utility::filePath(region_fasta_file, runtime_options_.workDirectory()) + ".fasta";
  kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 944950, 135, population_ptr_,
                                        get3D7Genome(), region_fasta_file);

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_941875";
  region_fasta_file = kgl::Utility::filePath(region_fasta_file, runtime_options_.workDirectory()) + ".fasta";
  kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 941875, 4293, population_ptr_,
                                        get3D7Genome(), region_fasta_file);

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_569342";
  region_fasta_file = kgl::Utility::filePath(region_fasta_file, runtime_options_.workDirectory()) + ".fasta";
  kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 569342, 7467, population_ptr_,
                                        get3D7Genome(), region_fasta_file);

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_584668";
  region_fasta_file = kgl::Utility::filePath(region_fasta_file, runtime_options_.workDirectory()) + ".fasta";
  kgl::GeneAnalysis::mutateGenomeRegion("malawi_fb_SRR609075", "Pf3D7_04_v3", 584668, 7280, population_ptr_,
                                        get3D7Genome(), region_fasta_file);


}


void kgl::PhylogeneticAnalysis::performMix() {

  // Index variants by ref/alt count.
  VariantClassifier classifier(unphased_population_ptr_);

  for (auto genome : classifier.getGenomes()) {

    const size_t minimum_offset_base_count = 25;
    const size_t maximum_offset_base_count = 150;
    kgd::MixtureDataObj mixture_data = classifier.convertToMixture(genome, minimum_offset_base_count,
                                                                   maximum_offset_base_count);
    kgd::ExecEnv::executeLib(mixture_data);

  }

}


void kgl::PhylogeneticAnalysis::performSNP() {


  // Split into country populations.
  std::string aux_file_path;
  if (not runtime_options_.getPropertiesAuxFile(aux_file_path)) {

    ExecEnv::log().critical("Aux Genome information file not found");

  }

  std::vector<CountryPair> country_pairs = GenomeAuxData::getCountries(aux_file_path, population_ptr_);
  GenomeAuxData aux_data;
  aux_data.readParseAuxData(aux_file_path);

  std::string DNA_mutation_file = kgl::Utility::filePath("DNAMutations.csv", runtime_options_.workDirectory());
  ApplicationAnalysis::outputDNAMutationCSV(DNA_mutation_file,
                                            PFATP4_CONTIG,
                                            PFATP4_GENE,
                                            PFATP4_SEQUENCE,
                                            get3D7Genome(),
                                            population_ptr_,
                                            aux_data);
  std::string amino_mutation_file = kgl::Utility::filePath("AminoMutations.csv", runtime_options_.workDirectory());
  ApplicationAnalysis::outputAminoMutationCSV(amino_mutation_file,
                                              PFATP4_CONTIG,
                                              PFATP4_GENE,
                                              PFATP4_SEQUENCE,
                                              get3D7Genome(),
                                              population_ptr_);


}


void kgl::PhylogeneticAnalysis::performUPGMA() {


  std::string newick_file = kgl::Utility::filePath("newick_VAR", runtime_options_.workDirectory()) + ".txt";
  std::string intron_file = kgl::Utility::filePath("intron_VAR", runtime_options_.workDirectory()) + ".csv";

  std::shared_ptr<const DNASequenceDistance> distance_metric_ptr(std::make_shared<const LevenshteinLocal>());

  UPGMAGeneFamilyTree<kgl::ReferenceGeneDistance>(newick_file,
                                                  intron_file,
                                                  distance_metric_ptr,
                                                  get3D7Genome(),
                                                  kgl::ReferenceGeneDistance::SYMBOLIC_VAR_FAMILY);

//  std::shared_ptr<const AminoSequenceDistance> distance_metric_ptr(std::make_shared<const Blosum80Global>());
//  kgl::UPGMAGenePhyloTree<kgl::UPGMAATP4Distance>(runtime_options_.workDirectory(),
//                                                  newick_file,
//                                                  distance_metric_ptr,
//                                                  population_ptr_,
//                                                  genome_db_ptr_,
//                                                  kgl::UPGMAProteinDistance::SYMBOLIC_VAR_FAMILY);

//  std::shared_ptr<const UnphasedPopulation> filtered_ptr = unphased_population_ptr->filterVariants(DPCountFilter(25));
//  kgl::UPGMAUnphasedTree<kgl::UPGMAUnphasedDistance>(newick_file, filtered_ptr);


}



void kgl::PhylogeneticAnalysis::performFineStructure() {


#define BASES_PER_CENTIMORGAN 15000.0

  // Split into country populations.
  std::string aux_file_path;
  if (not runtime_options_.getPropertiesAuxFile(aux_file_path)) {

    ExecEnv::log().critical("Aux Genome information file not found");

  }

  std::vector<CountryPair> country_pairs = GenomeAuxData::getCountries(aux_file_path, population_ptr_);

  for (auto country : country_pairs) {

    std::string file_name = "FS/" + country.first;
    std::string fine_structure_file = kgl::Utility::filePath(file_name, runtime_options_.workDirectory());
    FineStructureAnalysis::generateFiles(fine_structure_file, country.second, BASES_PER_CENTIMORGAN);

  }


}






void kgl::PhylogeneticAnalysis::performRNA() {


  RNAAnalysis rna_analysis;

  rna_analysis.getRNARegions("Pf3D7_04_v3",
                             956848,
                             135,
                             StrandSense::FORWARD,
                             "Pf3D7_04_v3",
                             (958066 - 1000),
                             9000,
                             StrandSense::FORWARD,
                             get3D7Genome());

  std::shared_ptr<const LocalDNASequenceCompare>  compare_metric_ptr(std::make_shared<const DNALocalAffineGap>());
  rna_analysis.compareRNARegion(0, 24, 1, compare_metric_ptr);

  rna_analysis.showResults(2);


}


