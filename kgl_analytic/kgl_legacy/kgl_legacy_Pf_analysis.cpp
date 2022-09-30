//
// Created by kellerberrin on 12/02/18.
//


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lagacy P.faciparum code. No longer used, waiting to be re-coded into analysis packages.


#include <kel_utility.h>
#include "kgl_sequence_distance.h"
#include "kgl_pfgenome_aux.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_analysis_gene_sequence.h"
#include "kgl_rna_search.h"
#include "kgl_legacy_Pf_analysis.h"

#include <functional>

namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Legacy P.Falciparum code. To be removed when the generalization to packages is complete.

/*
// Create a collection of all active organism genomes.
std::shared_ptr<kgl::GenomeCollection> genome_collection = GenomeCollection::createGenomeCollection(runtime_options);

// Create a phased population object.
std::shared_ptr<PhasedPopulation> population_ptr(std::make_shared<PhasedPopulation>("Organism"));

// Create an unphased population object.
std::shared_ptr<UnphasedPopulation> unphased_population_ptr(std::make_shared<UnphasedPopulation>("ParsedVCF"));

// Write Filtered Unphased Heterozygous Statistics
HeterozygousStatistics heterozygous_statistics;

RuntimeVCFFileMap vcf_map = runtime_options.getVCFFiles();
  // For all VCF files, read in the variants.
for (const auto &[vcf_ident, vcf_file] : vcf_map) {

  // Get VCF reference genome.
  std::shared_ptr<const GenomeReference> reference_genome_ptr = genome_collection->getGenomeGenealogyRecord(
  vcf_file.referenceGenome());

  // Filter and process Gatk variants.
  if (vcf_file.parserType() == VCFParserEnum::GatkMultiGenome) {

    // Read variants.
    std::shared_ptr<UnphasedPopulation> parsed_variants = VcfFactory::gatkMultiGenomeVCFVariants(
    reference_genome_ptr, vcf_file.fileName(), runtime_options.getEvidenceMap());

    // Basic statistics to output
    // unphased_population_ptr->popStatistics();
    // Filter unphased variants for minimum read statistics.
    std::shared_ptr<UnphasedPopulation> filtered_unphased_ptr = parsed_variants->filterVariants(
    AndFilter(DPCountFilter(30), RefAltCountFilter(30)));

    // Process Filtered Unphased Heterozygous Statistics
    if (not heterozygous_statistics.heterozygousStatistics(filtered_unphased_ptr)) {

      ExecEnv::log().error(
      "GeneExecEnv::executeApp(), Cannot generate heterozygous statistics for VCF file: {}",
      vcf_file.fileName());

    }

    // If the mixture file is defined and exists then read it and generate a population of clonal genomes.
    std::string mixture_file;
    if (runtime_options.getMixtureFile(mixture_file)) {

      // The mixture file (Pf3k only) indicates the Complexity Of Infection (COI) of VCF samples, this function includes only clonal infections.
      std::shared_ptr<UnphasedPopulation> clonal_unphased = GenomePhasing::filterClonal(mixture_file, filtered_unphased_ptr);

    }

    // Get the VCF ploidy (need not be the organism ploidy). Only applies to Pf3k, removed from VCF specification.
    size_t ploidy = 2;
    // Phase the homozygous and heterozygous variants into a haploid population.
    GenomePhasing::haploidPhasing(ploidy, filtered_unphased_ptr, reference_genome_ptr, population_ptr);


    std::pair<size_t, size_t> valid_count = filtered_unphased_ptr->validate(reference_genome_ptr);
    ExecEnv::log().info("Population: {}, Total Variants: {}, Validated Variants: {}", filtered_unphased_ptr->populationId(), valid_count.first, valid_count.second);

    unphased_population_ptr->mergePopulation(filtered_unphased_ptr);

  } else if (vcf_file.parserType() == VCFParserEnum::GRChNoGenome) {

    ContigAliasMap contig_alias_map = runtime_options.getContigAlias();
    // Read variants.
    std::shared_ptr<UnphasedGenome> parsed_variants = VcfFactory::GRChNoGenomeVCFVariants(reference_genome_ptr,
                                                                                          vcf_file.fileName(),
                                                                                          contig_alias_map,
                                                                                          runtime_options.getEvidenceMap());

    std::pair<size_t, size_t> valid_count = parsed_variants->validate(reference_genome_ptr);
    ExecEnv::log().info("Genome: {}, Total Variants: {}, Validated Variants: {}", parsed_variants->genomeId(), valid_count.first, valid_count.second);

  }

}


// Write the unphased hetero/homozygous statistics.
std::string heterozygous_file = Utility::filePath("HeterozygousAnalysis", args.workDirectory) + ".csv";
if (not heterozygous_statistics.writeHeterozygousStatistics(heterozygous_file, ',')) {

  ExecEnv::log().error("GeneExecEnv::executeApp(), Cannot write heterozygous statistics to file: {}", heterozygous_file);

}

// Analyze the data.
kgl::PhylogeneticAnalysis(runtime_options, genome_collection, unphased_population_ptr, population_ptr).performAnalysis("INTERVAL");

}


*/


const kgl::GenomeId_t analysis_genome = "Pf3D7_47";


void kgl::PhylogeneticAnalysis::performAnalysis(const std::string& analysis_type) {

std::vector<std::string> analysis_ids = Utility::tokenizer(analysis_type, ANALYSIS_SEPARATORS_);

if (analysis_ids.empty()) {

  ExecEnv::log().info("No Analysis functions specified, valid analysis below:");

  for (auto const& analysis: analysis_map_) {

    ExecEnv::log().info("Valid analysis type: '{}'", analysis.first);

  }

} else {

  for (auto const& analysis : analysis_ids) {

    dispatchAnalysis(analysis);

  }

}

}

void kgl::PhylogeneticAnalysis::dispatchAnalysis(const std::string& analysis_type) {

auto result = analysis_map_.find(analysis_type);

if (result == analysis_map_.end()) {

  ExecEnv::log().error("Invalid Analysis Type: '{}'. Must be one of the following (case insensitive):", analysis_type);

  for (auto analysis: analysis_map_) {

    ExecEnv::log().info("Valid Analysis Type: '{}'", analysis.first);

  }

  return;

}

ExecEnv::log().info("**********Performing Analysis: {} **************", result->first);

// Call the analysis.
std::invoke(result->second, this);

ExecEnv::log().info("**********Completed Analysis: {} **************", result->first);

}



void kgl::PhylogeneticAnalysis::performSequence() {


  std::shared_ptr<const DNASequenceDistance> dna_distance_metric(std::make_shared<const LevenshteinGlobal>());

  std::string coding_file = Utility::filePath("All_DNA_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  if (not kgl::GenomicMutation::outputDNASequenceCSV(coding_file,
                                                     SequenceAnalysisType::DNA,
                                                     dna_distance_metric,
                                                     genome_collection_ptr_->getGenome(analysis_genome),
                                                     population_ptr_)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

  }

  coding_file = Utility::filePath("All_Variant_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  if (not kgl::GenomicMutation::outputDNASequenceCSV(coding_file, SequenceAnalysisType::VARIANT, dna_distance_metric,
                                                     genome_collection_ptr_->getGenome(analysis_genome),
                                                     population_ptr_)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

  }

  coding_file = Utility::filePath("All_SNP_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  if (not kgl::GenomicMutation::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SNP, dna_distance_metric,
                                                     genome_collection_ptr_->getGenome(analysis_genome),
                                                     population_ptr_)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

  }

  coding_file = Utility::filePath("All_SIZE_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  if (not kgl::GenomicMutation::outputDNASequenceCSV(coding_file, SequenceAnalysisType::SIZE, dna_distance_metric,
                                                     genome_collection_ptr_->getGenome(analysis_genome),
                                                     population_ptr_)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

  }

  std::shared_ptr<const AminoSequenceDistance> amino_distance_metric(std::make_shared<const LevenshteinGlobal>());

  coding_file = Utility::filePath("All_Amino_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
  if (not kgl::GenomicMutation::outputAminoSequenceCSV(coding_file, amino_distance_metric,
                                                       genome_collection_ptr_->getGenome(analysis_genome),
                                                       population_ptr_)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

  }

  // Split into country populations.
  std::string aux_file_path;

  std::vector<PfCountryPair> country_pairs = PfGenomeAuxData::getCountries(aux_file_path, population_ptr_);

  for (auto country : country_pairs) {

    coding_file = Utility::filePath(country.first + "Variant_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    if (not kgl::GenomicMutation::outputDNASequenceCSV(coding_file,
                                                       SequenceAnalysisType::VARIANT, dna_distance_metric,
                                                       genome_collection_ptr_->getGenome(analysis_genome),
                                                       country.second)) {

      ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

    }

    coding_file = Utility::filePath(country.first + "SNP_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    if (not kgl::GenomicMutation::outputDNASequenceCSV(coding_file,
                                                       SequenceAnalysisType::SNP, dna_distance_metric,
                                                       genome_collection_ptr_->getGenome(analysis_genome),
                                                       country.second)) {

      ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

    }

    coding_file = Utility::filePath(country.first + "SIZE_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    if (not kgl::GenomicMutation::outputDNASequenceCSV(coding_file,
                                                       SequenceAnalysisType::SIZE, dna_distance_metric,
                                                       genome_collection_ptr_->getGenome(analysis_genome),
                                                       country.second)) {

      ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

    }

    coding_file = Utility::filePath(country.first + "DNA_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    if (kgl::GenomicMutation::outputDNASequenceCSV(coding_file,
                                                   SequenceAnalysisType::DNA, dna_distance_metric,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   country.second)) {

      ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

    }

    coding_file = Utility::filePath(country.first + "Amino_CodingAnalysis", runtime_options_.workDirectory()) + ".csv";
    if (not kgl::GenomicMutation::outputAminoSequenceCSV(coding_file,
                                                         amino_distance_metric,
                                                         genome_collection_ptr_->getGenome(analysis_genome),
                                                         country.second)) {

      ExecEnv::log().error("PhylogeneticAnalysis::performSequence(), outputDNASequenceCSV fails");

    }

  }

}



void kgl::PhylogeneticAnalysis::performGene() {

std::string fasta_file = Utility::filePath(ACTIVE_SEQUENCE, runtime_options_.workDirectory()) + ".fasta";

if (not kgl::GenomicSequence::translateContig("NC_045512_2", "NC_045512.2", genome_collection_ptr_, fasta_file)) {

  ExecEnv::log().error("GenomicSequence::translateContig(); analysis fails");

}

}


void kgl::PhylogeneticAnalysis::performMix() {


}


void kgl::PhylogeneticAnalysis::performSNP() {


// Split into country populations.
std::string aux_file_path;

std::vector<PfCountryPair> country_pairs = PfGenomeAuxData::getCountries(aux_file_path, population_ptr_);
PfGenomeAuxData aux_data;
if (not aux_data.readParseAuxData(aux_file_path)) {

  ExecEnv::log().critical("performSNP; Cannot read/parse aux file: {}", aux_file_path);

}


std::string DNA_mutation_file = Utility::filePath("DNAMutations.csv", runtime_options_.workDirectory());
if (not GenomicMutation::outputDNAMutationCSV(DNA_mutation_file,
                                              PFATP4_CONTIG,
                                              PFATP4_GENE,
                                              PFATP4_SEQUENCE,
                                              genome_collection_ptr_->getGenome(analysis_genome),
                                              population_ptr_,
                                              aux_data)) {

  ExecEnv::log().error("PhylogeneticAnalysis::performSNP(), outputDNAMutationCSV fails");

}
std::string amino_mutation_file = Utility::filePath("AminoMutations.csv", runtime_options_.workDirectory());
if (not GenomicMutation::outputAminoMutationCSV(amino_mutation_file,
                                                PFATP4_CONTIG,
                                                PFATP4_GENE,
                                                PFATP4_SEQUENCE,
                                                genome_collection_ptr_->getGenome(analysis_genome),
                                                population_ptr_)) {

  ExecEnv::log().error("PhylogeneticAnalysis::performSNP(), outputDNAMutationCSV fails");

}


}






void kgl::PhylogeneticAnalysis::performRNA() {


RNAAnalysis rna_analysis;

if (not rna_analysis.getRNARegions("Pf3D7_04_v3",
                           956848,
                           135,
                           StrandSense::FORWARD,
                           "Pf3D7_04_v3",
                           (958066 - 1000),
                           9000,
                           StrandSense::FORWARD,
                                   genome_collection_ptr_->getGenome(analysis_genome))) {

  ExecEnv::log().error("Problem performing RNA analysis");

}

std::shared_ptr<const LocalDNASequenceCompare>  compare_metric_ptr(std::make_shared<const DNALocalAffineGap>());

if (not rna_analysis.compareRNARegion(0, 24, 1, compare_metric_ptr)) {

  ExecEnv::log().error("Problem performing RNA analysis");

}

rna_analysis.showResults(2);


}


