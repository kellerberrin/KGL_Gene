//
// Created by kellerberrin on 16/01/18.
//

#include "kgl_sequence_offset.h"
#include "kgl_phylogenetic_gene.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_sequence_complexity.h"

#include <memory>


namespace kgl = kellerberrin::genome;


bool kgl::GeneAnalysis::translateGene( const GenomeId_t& genome_id, 
                                       const FeatureIdent_t& gene_id,
                                       std::shared_ptr<const GenomeCollection>& genomes, 
                                       const std::string& fasta_file_name) {

  // Get the genome object
  std::shared_ptr<const GenomeDatabase> genome_db_ptr;
  if (not genomes->getGenome(genome_id, genome_db_ptr)) {

    ExecEnv::log().warn("Could not find Genome: {} in genome collection", genome_id);
    return false;
    
  }

  // Look through the contigs for the gene.
  bool found = false;
  std::vector<std::shared_ptr<const Feature>> feature_ptr_vec;
  for (const auto& contig : genome_db_ptr->getMap() ) {

    if (contig.second->findFeatureId(gene_id, feature_ptr_vec)) {

      found = true;
      break;   

    }

  }

  if (not found or feature_ptr_vec.empty()) {

    ExecEnv::log().warn("Did not find Gene: {} in Genome: {}", gene_id, genome_id);
    return false;

  }

  std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>> amino_fasta_vector;
  for (const auto& feature_ptr : feature_ptr_vec) {

    
    if (feature_ptr->isGene()) {

      std::shared_ptr<const kgl::CodingSequenceArray> coding_sequence_array_ptr;
      coding_sequence_array_ptr = GeneFeature::getCodingSequences(std::static_pointer_cast<const GeneFeature>(feature_ptr)); 

      for (const auto& coding_sequence : coding_sequence_array_ptr->getMap()) {

        std::shared_ptr<DNA5SequenceCoding> coding_dna_ptr;
        if (coding_sequence.second->contig()->getDNA5SequenceCoding(coding_sequence.second, coding_dna_ptr)) {

          std::shared_ptr<AminoSequence> peptide_ptr = coding_sequence.second->contig()->getAminoSequence(coding_dna_ptr);
          std::string reference_name = gene_id;
          std::pair<std::string, std::shared_ptr<AminoSequence>> fasta_entry(reference_name, peptide_ptr);
          amino_fasta_vector.push_back(fasta_entry);

        } else {

          ExecEnv::log().error("Cannot generate a DNA coding sequence for Gene: {} in Genome: {}", gene_id, genome_id);
          return false;

        }

      }

    }

  }

  // Sequences are presented as a pair of a sequence name and an amino sequences.
  if (not ApplicationAnalysis::writeMutantProteins(fasta_file_name, amino_fasta_vector)) {

    ExecEnv::log().warn("Problem writing peptide fasta for Gene: {} in file: {}", gene_id, fasta_file_name);
    return false;

  }

  return true;

}



bool kgl::GeneAnalysis::mutateGene(const ContigId_t& contig,
                                   const FeatureIdent_t& gene,
                                   const FeatureIdent_t& sequence,
                                   std::shared_ptr<const PhasedPopulation> population_ptr,
                                   std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                   const std::string& fasta_filename) {



  GeneSummaryMap gene_summary_map;

  for (auto genome : population_ptr->getMap()) {

    if (not mutateGenomeGene(contig, gene, sequence, genome.second, genome_db_ptr, gene_summary_map)) {

      return false;

    }

  }

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db_ptr->getContigSequence(contig, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig);
    return false;

  }

  for (auto summary : gene_summary_map) {

    if (contig_ptr->verifyProteinSequence(summary.second.sequence_mutant)) {

      ExecEnv::log().info("Multiple comparison Genome: {}, protein score: {}, 5Prime score: {}, 3Prime score: {}",
                          summary.first, summary.second.sequence_distance, summary.second.prime5_distance, summary.second.prime3_distance);

    } else {

      ExecEnv::log().info("Frame shift mutation Genome: {}, protein score: {}, 5Prime score: {}, 3Prime score: {}",
                          summary.first, summary.second.sequence_distance, summary.second.prime5_distance, summary.second.prime3_distance);


    }

  }

  // Write all the amino sequences as a fasta file.
  std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>> amino_fasta_vector;
  for (auto summary : gene_summary_map) {

    std::string reference_name = summary.first + "_Reference_" + sequence;
    std::pair<std::string, std::shared_ptr<AminoSequence>> fasta_entry(reference_name, summary.second.sequence_ptr);
    amino_fasta_vector.push_back(fasta_entry);

    size_t mutant_count = 0;

    ++mutant_count;
    std::stringstream ss;
    ss << summary.first << "_Mutant" << mutant_count << "_" << sequence;
    std::pair<std::string, std::shared_ptr<AminoSequence>> mutant_fasta_entry(ss.str(), summary.second.sequence_mutant);
    amino_fasta_vector.push_back(mutant_fasta_entry);

  }

  // Sequences are presented as a pair of a sequence name and an amino sequences.
  if (not ApplicationAnalysis::writeMutantProteins(fasta_filename, amino_fasta_vector)) {

    ExecEnv::log().warn("mutantProtein(), Problem writing protein fasta file: {}", fasta_filename);

  }


  return true;

}



bool kgl::GeneAnalysis::mutateGenomeGene(const ContigId_t& contig,
                                         const FeatureIdent_t& gene,
                                         const FeatureIdent_t& sequence,
                                         std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                         GeneSummaryMap& gene_summary_map) {


  GeneSummary gene_summary;
  std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric(std::make_shared<const LevenshteinGlobal>());
  std::shared_ptr<const GlobalAminoSequenceDistance> amino_distance_metric(std::make_shared<const LevenshteinGlobal>());
  gene_summary.genome = genome_variant_ptr->genomeId();
  if (ApplicationAnalysis::compare5Prime(contig,
                                         gene,
                                         sequence,
                                         PRIME_REGION_SIZE,
                                         genome_db_ptr,
                                         genome_variant_ptr,
                                         gene_summary.prime5_reference,
                                         gene_summary.prime5_mutant)) {

    gene_summary.prime5_distance = dna_distance_metric->distance(gene_summary.prime5_reference, gene_summary.prime5_mutant);
    ExecEnv::log().info("5PRIME Genome: {}, Contig: {}, Gene: {}, Sequence: {} Levenshtein: {}, comparison:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.prime5_distance);

  } else {

    ExecEnv::log().error("MutateGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                         genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }


  if (genome_variant_ptr->mutantProteins(contig,
                                         ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                         gene,
                                         sequence,
                                         genome_db_ptr,
                                         gene_summary.variant_map,
                                         gene_summary.sequence_ptr,
                                         gene_summary.sequence_mutant)) {

    std::stringstream ss;
    ss << gene_summary.variant_map;

    ExecEnv::log().info("Variants used to mutate Genome: {}, Contig: {}, Gene: {} Sequence: {}:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, ss.str());

    gene_summary.sequence_distance = amino_distance_metric->distance( gene_summary.sequence_ptr, gene_summary.sequence_mutant);
    ExecEnv::log().info("Genome: {}, Contig: {}, Gene: {}, Sequence: {} Levenshtein: {}, comparison:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.sequence_distance);

  } else {

    ExecEnv::log().error("MutateGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }

  if (ApplicationAnalysis::compare3Prime(contig,
                                         gene,
                                         sequence,
                                         PRIME_REGION_SIZE,
                                         genome_db_ptr,
                                         genome_variant_ptr,
                                         gene_summary.prime3_reference,
                                         gene_summary.prime3_mutant)) {

    gene_summary.prime3_distance = dna_distance_metric->distance(gene_summary.prime3_reference, gene_summary.prime3_mutant);
    ExecEnv::log().info("3PRIME Genome: {}, Contig: {}, Gene: {}, Sequence: {} Levenshtein: {}, comparison:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.prime3_distance);

  } else {

    ExecEnv::log().error("MutateGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                         genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }

  auto result = gene_summary_map.insert(std::pair<GenomeId_t, GeneSummary>(gene_summary.genome, gene_summary));

  if (not result.second) {

    ExecEnv::log().error("MutateGenomeGene(), Duplicate sequence; Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                         genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }


  return true;

}


bool kgl::GeneAnalysis::mutateAllRegions(const std::string& file_name,
                                         ContigSize_t region_size,
                                         std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric,
                                         std::shared_ptr<const PhasedPopulation> pop_variant_ptr,
                                         std::shared_ptr<const GenomeDatabase> genome_db_ptr) {

  const char CSV_delimiter = ',';
  // open the file.
  std::fstream out_file(file_name, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  out_file << outputRegionHeader(CSV_delimiter) << '\n';

  for( auto genome_variant : pop_variant_ptr->getMap()) {

    ExecEnv::log().info("outputSequenceCSV(), Processing genome: {}", genome_variant.first);
    ContigSize_t sequence_count = 0;

    for (auto contig : genome_db_ptr->getMap()) {

      for (ContigOffset_t offset = 0; offset < (contig.second->sequence().length() - region_size); offset += region_size) {

        out_file << outputGenomeRegion(CSV_delimiter, dna_distance_metric, contig.first, offset, region_size, genome_variant.second, genome_db_ptr) << '\n';
        ++sequence_count;

      }

    }

    ExecEnv::log().info("outputSequenceCSV(), Genome: {} mutated: {} sequences.", genome_variant.first, sequence_count);

  }

  return out_file.good();

}


std::string kgl::GeneAnalysis::outputRegionHeader(char delimiter) {

  std::stringstream ss;

  ss << "Genome" << delimiter;
  ss << "Contig" << delimiter;
  ss << "ContigSize" << delimiter;
  ss << "ContigOffset" << delimiter;
  ss << "RegionSize" << delimiter;
  ss << "Score" << delimiter;
  ss << "CpG" << delimiter;
  ss << "A_prop" << delimiter;
  ss << "C_prop" << delimiter;
  ss << "G_prop" << delimiter;
  ss << "T_prop";

  return ss.str();

}



std::string kgl::GeneAnalysis::outputGenomeRegion(char delimiter,
                                                  std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric,
                                                  const ContigId_t& contig_id,
                                                  const ContigOffset_t offset,
                                                  const ContigSize_t region_size,
                                                  std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr) {

  std::stringstream ss;

  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db_ptr->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().error("outputGenomeGene(), Unexpected could not find Contig: {}", contig_id);
    return "<error>";

  }
  if (offset + region_size >= contig_ptr->sequence().length()) {

    ExecEnv::log().error("outputGenomeGene(), Offset: {} + Region Size: {} exceed the Contig: {} size: {}",
                         offset, region_size, contig_id, contig_ptr->sequence().length());
    return "<error>";

  }


  if (region_size == 0) {

    ss << genome_variant_ptr->genomeId() << delimiter;
    ss << contig_id << delimiter;
    ss << contig_ptr->sequence().length() << delimiter;
    ss << offset << delimiter;
    ss << region_size << delimiter;
    ss << 0 << delimiter;
    ss << 0.0 << delimiter;
    ss << 0.0 << delimiter;
    ss << 0.0 << delimiter;
    ss << 0.0 << delimiter;
    ss << 0.0;

    return ss.str();

  }

  std::shared_ptr<DNA5SequenceLinear> mutant_sequence;
  std::shared_ptr<DNA5SequenceLinear> reference_sequence;
  OffsetVariantMap variant_map;
  if (genome_variant_ptr->mutantRegion(contig_id,
                                       ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                       offset,
                                       region_size,
                                       genome_db_ptr,
                                       variant_map,
                                       reference_sequence,
                                       mutant_sequence)) {


    double distance = static_cast<double>(dna_distance_metric->distance(reference_sequence, mutant_sequence));
    ss << genome_variant_ptr->genomeId() << delimiter;
    ss << contig_id << delimiter;
    ss << contig_ptr->sequence().length() << delimiter;
    ss << offset << delimiter;
    ss << region_size << delimiter;
    ss << distance << delimiter;
    ss << SequenceComplexity::relativeCpGIslands(reference_sequence) << delimiter;
    double A_prop;
    double C_prop;
    double G_prop;
    double T_prop;
    SequenceComplexity::proportionNucleotides(reference_sequence, A_prop, C_prop, G_prop, T_prop);
    ss << A_prop << delimiter;
    ss << C_prop << delimiter;
    ss << G_prop << delimiter;
    ss << T_prop;

  } else {

    ExecEnv::log().error("outputGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Offset: {}, Size: {}",
                         genome_variant_ptr->genomeId(), contig_id, offset, region_size);
    return "<error>";

  }

  return ss.str();

}



bool kgl::GeneAnalysis::mutateGenomeRegion(const GenomeId_t& genome,
                                           const ContigId_t& contig,
                                           ContigOffset_t offset,
                                           ContigSize_t region_size,
                                           std::shared_ptr<const PhasedPopulation> population_ptr,
                                           std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                           const std::string& fasta_file) {


  std::shared_ptr<const GenomeVariant> genome_variant_ptr;
  if (population_ptr->getGenomeVariant(genome, genome_variant_ptr)) {

    return mutateGenomeRegion(contig, offset, region_size, genome_variant_ptr, genome_db_ptr, fasta_file);

  } else {

    ExecEnv::log().error("mutateGenomeRegion(), Genome: {} not found", genome);
    return false;

  }

}


bool kgl::GeneAnalysis::mutateGenomeRegion(const ContigId_t& contig,
                                           const ContigOffset_t offset,
                                           const ContigSize_t region_size,
                                           std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                                           std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                           const std::string& fasta_file) {

  std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric(std::make_shared<const LevenshteinGlobal>());
  std::vector<std::pair<std::string, std::shared_ptr<VirtualSequence>>> dna_seq_vector;
  std::shared_ptr<DNA5SequenceLinear> mutant_sequence;
  std::shared_ptr<DNA5SequenceLinear> reference_sequence;
  OffsetVariantMap variant_map;
  if (genome_variant_ptr->mutantRegion(contig,
                                       ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                       offset,
                                       region_size,
                                       genome_db_ptr,
                                       variant_map,
                                       reference_sequence,
                                       mutant_sequence)) {

    std::stringstream ss;
    ss << variant_map;

    ExecEnv::log().info("Variants used to mutate Genome: {}, Contig: {}, Offset: {} Size: {}:\n{}",
                        genome_variant_ptr->genomeId(), contig, offset, region_size, ss.str());

    std::stringstream ref_fasta_ss;
    ref_fasta_ss << "reference_" << contig << "_+_" << offset << "_" << region_size;
    std::pair<std::string, std::shared_ptr<VirtualSequence>> fasta_reference(ref_fasta_ss.str(), reference_sequence);
    dna_seq_vector.push_back(fasta_reference);

    std::shared_ptr<DNA5SequenceCoding> ref_reverse_complement_ptr = SequenceOffset::codingSequence(reference_sequence, StrandSense::REVERSE);

    ref_fasta_ss.str("");
    ref_fasta_ss << "reference_" << contig << "_-_" << offset << "_" << region_size;
    std::pair<std::string, std::shared_ptr<VirtualSequence>> rev_fasta_reference(ref_fasta_ss.str(), ref_reverse_complement_ptr);
    dna_seq_vector.push_back(rev_fasta_reference);


    size_t count = 0;

    CompareDistance_t score = dna_distance_metric->distance(reference_sequence, mutant_sequence);
    ExecEnv::log().info("Genome: {}, Contig: {}, Offset: {} Size: {} score: {}",
                        genome_variant_ptr->genomeId(), contig, offset, region_size, score);

    count++;
    std::stringstream mutant_fasta_ss;
    mutant_fasta_ss << "mutant_" << count << "_" << contig << "_" << offset << "_" << region_size;
    std::pair<std::string, std::shared_ptr<DNA5SequenceLinear>> fasta_mutant(mutant_fasta_ss.str(), mutant_sequence);
    dna_seq_vector.push_back(fasta_mutant);

    if (not ParseGffFasta().writeFastaFile(fasta_file, dna_seq_vector)) {

      ExecEnv::log().error("MutateGenomeGene(), Problem writing to fasta file: {}", fasta_file);

    }

  } else {

    ExecEnv::log().error("MutateGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Offset: {}, Size: {}",
                         genome_variant_ptr->genomeId(), contig, offset, region_size);
    return false;

  }

  return true;

}

