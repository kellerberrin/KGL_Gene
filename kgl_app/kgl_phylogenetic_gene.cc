//
// Created by kellerberrin on 16/01/18.
//

#include "kgl_sequence_offset.h"
#include "kgl_phylogenetic_gene.h"
#include "kgl_phylogenetic_analysis.h"


namespace kgl = kellerberrin::genome;


bool kgl::GeneAnalysis::mutateGene(const ContigId_t& contig,
                                   const FeatureIdent_t& gene,
                                   const FeatureIdent_t& sequence,
                                   std::shared_ptr<const PopulationVariant> population_ptr,
                                   std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                   const std::string& fasta_filename) {



  GeneSummaryMap gene_summary_map;

  for (auto genome : population_ptr->getMap()) {

    if (not mutateGenomeGene(contig, gene, sequence, genome.second, genome_db_ptr, gene_summary_map)) {

      return false;

    }

  }

  std::vector<std::shared_ptr<const DNA5SequenceCoding>> prime5_vector;
  for (auto summary : gene_summary_map) {

    for (auto mutant : summary.second.prime5_mutant_vec) {

      prime5_vector.push_back(mutant);

    }

  }

  std::string compare = DNA5SequenceCoding::multipleCompare(prime5_vector);
  ExecEnv::log().info("Comparison of all 5 prime regions for Contig: {}, Gene: {}, Sequence: {}; \n{}",
                      contig, gene, sequence, compare);


  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db_ptr->getContigSequence(contig, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig);
    return false;

  }

  std::vector<std::shared_ptr<const AminoSequence>> amino_vector;
  for (auto summary : gene_summary_map) {

    for (auto mutant : summary.second.sequence_mutant_vec) {

      if (contig_ptr->verifyProteinSequence(mutant)) {

        amino_vector.push_back(mutant);
        ExecEnv::log().info("Multiple comparison Genome: {}, protein score: {}, 5Prime score: {}, 3Prime score: {}",
                            summary.first, summary.second.sequence_score, summary.second.prime5_score, summary.second.prime3_score);

      } else {

        ExecEnv::log().info("Frame shift mutation Genome: {}, protein score: {}, 5Prime score: {}, 3Prime score: {}",
                            summary.first, summary.second.sequence_score, summary.second.prime5_score, summary.second.prime3_score);


      }

    }

  }



  // Write all the amino sequences as a fasta file.
  std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>> amino_fasta_vector;
  for (auto summary : gene_summary_map) {

    std::string reference_name = summary.first + "_Reference_" + sequence;
    std::pair<std::string, std::shared_ptr<AminoSequence>> fasta_entry(reference_name, summary.second.sequence_ptr);
    amino_fasta_vector.push_back(fasta_entry);

    size_t mutant_count = 0;
    for (auto mutant : summary.second.sequence_mutant_vec) {

      ++mutant_count;
      std::stringstream ss;
      ss << summary.first << "_Mutant" << mutant_count << "_" << sequence;
      std::pair<std::string, std::shared_ptr<AminoSequence>> mutant_fasta_entry(ss.str(), mutant);
      amino_fasta_vector.push_back(mutant_fasta_entry);

    }

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
  gene_summary.genome = genome_variant_ptr->genomeId();
  if (ApplicationAnalysis::compare5Prime(contig,
                                         gene,
                                         sequence,
                                         PRIME_REGION_SIZE,
                                         genome_db_ptr,
                                         genome_variant_ptr,
                                         gene_summary.prime5_reference,
                                         gene_summary.prime5_mutant_vec)) {

    for (auto mutant : gene_summary.prime5_mutant_vec) {

      std::string comparison = gene_summary.prime5_reference->compareDNA5Coding(mutant, gene_summary.prime5_score);
      ExecEnv::log().info("5PRIME Genome: {}, Contig: {}, Gene: {}, Sequence: {} score: {}, comparison:\n{}",
                          genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.prime5_score, comparison);

    }

  } else {

    ExecEnv::log().error("MutateGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                         genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }


  if (genome_variant_ptr->mutantProteins(contig,
                                         gene,
                                         sequence,
                                         genome_db_ptr,
                                         gene_summary.variant_map,
                                         gene_summary.sequence_ptr,
                                         gene_summary.sequence_mutant_vec)) {

    std::stringstream ss;
    ss << gene_summary.variant_map;

    ExecEnv::log().info("Variants used to mutate Genome: {}, Contig: {}, Gene: {} Sequence: {}:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, ss.str());

    for (auto mutant : gene_summary.sequence_mutant_vec) {

      std::string comparison = gene_summary.sequence_ptr->compareAminoSequences(mutant, gene_summary.sequence_score);
      ExecEnv::log().info("Genome: {}, Contig: {}, Gene: {}, Sequence: {} score: {}, comparison:\n{}",
                          genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.sequence_score, comparison);

    }

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
                                         gene_summary.prime3_mutant_vec)) {

    for (auto mutant : gene_summary.prime3_mutant_vec) {

      std::string comparison = gene_summary.prime3_reference->compareDNA5Coding(mutant, gene_summary.prime3_score);
      ExecEnv::log().info("3PRIME Genome: {}, Contig: {}, Gene: {}, Sequence: {} score: {}, comparison:\n{}",
                          genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.prime3_score, comparison);

    }

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
                                         std::shared_ptr<const PopulationVariant> pop_variant_ptr,
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

        out_file << outputGenomeRegion(CSV_delimiter, contig.first, offset, region_size, genome_variant.second, genome_db_ptr) << '\n';
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
  ss << "RegionGC" << delimiter;
  ss << "ContigOffset" << delimiter;
  ss << "RegionSize" << delimiter;
  ss << "Score";

  return ss.str();

}


std::string kgl::GeneAnalysis::outputGenomeRegion(char delimiter,
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

  double proportion_GC = 0;

  if (region_size == 0) {

    ss << genome_variant_ptr->genomeId() << delimiter;
    ss << contig_id << delimiter;
    ss << contig_ptr->sequence().length() << delimiter;
    ss << proportion_GC << delimiter;
    ss << offset << delimiter;
    ss << region_size << delimiter;
    ss << 0;

    return ss.str();

  }

  std::vector<std::shared_ptr<DNA5SequenceLinear>> mutant_sequence_vector;
  std::shared_ptr<DNA5SequenceLinear> reference_sequence;
  OffsetVariantMap variant_map;
  if (genome_variant_ptr->mutantRegion(contig_id,
                                       offset,
                                       region_size,
                                       genome_db_ptr,
                                       variant_map,
                                       reference_sequence,
                                       mutant_sequence_vector)) {

    double average_score = 0;
    for (auto mutant : mutant_sequence_vector) {

      CompareScore_t score;
      std::string comparison = reference_sequence->compareDNA5Sequences(mutant, score);
      average_score += static_cast<double>(score);

    }

    proportion_GC = static_cast<double>(reference_sequence->countGC()) / static_cast<double>(region_size);
    average_score = average_score / static_cast<double>(mutant_sequence_vector.size());

    ss << genome_variant_ptr->genomeId() << delimiter;
    ss << contig_id << delimiter;
    ss << contig_ptr->sequence().length() << delimiter;
    ss << proportion_GC << delimiter;
    ss << offset << delimiter;
    ss << region_size << delimiter;
    ss << average_score;

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
                                           std::shared_ptr<const PopulationVariant> population_ptr,
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

  const size_t MAX_SEQUENCE_DISPLAY_SIZE = 25000; // largest sequence to display
  std::vector<std::pair<std::string, std::shared_ptr<VirtualSequence>>> dna_seq_vector;
  std::vector<std::shared_ptr<DNA5SequenceLinear>> mutant_sequence_vector;
  std::shared_ptr<DNA5SequenceLinear> reference_sequence;
  OffsetVariantMap variant_map;
  if (genome_variant_ptr->mutantRegion(contig,
                                       offset,
                                       region_size,
                                       genome_db_ptr,
                                       variant_map,
                                       reference_sequence,
                                       mutant_sequence_vector)) {

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
    for (auto mutant : mutant_sequence_vector) {

      if (reference_sequence->length() < MAX_SEQUENCE_DISPLAY_SIZE) {
        CompareScore_t score;
        std::string comparison = reference_sequence->compareDNA5Sequences(mutant, score);
        double normed_score = static_cast<double>(score) * 1000.0 / static_cast<double>(mutant->length());
        ExecEnv::log().info("Genome: {}, Contig: {}, Offset: {} Size: {} score: {} score/1000: {} comparison:\n{}",
                            genome_variant_ptr->genomeId(), contig, offset, region_size, score, normed_score, comparison);
      } else {

        CompareScore_t score = reference_sequence->compareLevenshtein(mutant);
        double normed_score = static_cast<double>(score) * 1000.0 / static_cast<double>(mutant->length());
        ExecEnv::log().info("Genome: {}, Contig: {}, Offset: {} Size: {} score: {} score/1000: {} \n ....Sequence too large to display (Levenshtein Sequence Match)....",
                            genome_variant_ptr->genomeId(), contig, offset, region_size, score, normed_score);

      }

      count++;
      std::stringstream mutant_fasta_ss;
      mutant_fasta_ss << "mutant_" << count << "_" << contig << "_" << offset << "_" << region_size;
      std::pair<std::string, std::shared_ptr<DNA5SequenceLinear>> fasta_mutant(mutant_fasta_ss.str(), mutant);
      dna_seq_vector.push_back(fasta_mutant);

    }

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

