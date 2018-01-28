//
// Created by kellerberrin on 16/01/18.
//

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

//  std::string amino_compare = AminoSequence::multipleCompare(amino_vector);
//  ExecEnv::log().info("Comparison of all Amino sequences for Contig: {}, Gene: {}, Sequence: {}; \n{}",
//                      contig, gene, sequence, amino_compare);


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



bool kgl::GeneAnalysis::mutateRegion(const ContigId_t& contig,
                                     ContigOffset_t offset,
                                     ContigSize_t region_size,
                                     std::shared_ptr<const PopulationVariant> population_ptr,
                                     std::shared_ptr<const GenomeDatabase> genome_db_ptr) {


  for (auto genome : population_ptr->getMap()) {

    if (not mutateGenomeRegion(contig, offset, region_size, genome.second, genome_db_ptr)) {

      return false;

    }

  }

  return true;

}

bool kgl::GeneAnalysis::mutateGenomeRegion(const GenomeId_t& genome,
                                           const ContigId_t& contig,
                                           ContigOffset_t offset,
                                           ContigSize_t region_size,
                                           std::shared_ptr<const PopulationVariant> population_ptr,
                                           std::shared_ptr<const GenomeDatabase> genome_db_ptr) {


  std::shared_ptr<const GenomeVariant> genome_variant_ptr;
  if (population_ptr->getGenomeVariant(genome, genome_variant_ptr)) {

    return mutateGenomeRegion(contig, offset, region_size, genome_variant_ptr, genome_db_ptr);

  } else {

    ExecEnv::log().error("mutateGenomeRegion(), Genome: {} not found", genome);
    return false;

  }

}


bool kgl::GeneAnalysis::mutateGenomeRegion(const ContigId_t& contig,
                                           const ContigOffset_t offset,
                                           const ContigSize_t region_size,
                                           std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                                           std::shared_ptr<const GenomeDatabase> genome_db_ptr) {




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


    for (auto mutant : mutant_sequence_vector) {

      CompareScore_t score;
      std::string comparison = reference_sequence->compareDNA5Sequences(mutant, score);
      ExecEnv::log().info("Genome: {}, Contig: {}, Offset: {} Size: {} score: {}, comparison:\n{}",
                          genome_variant_ptr->genomeId(), contig, offset, region_size, score, comparison);

    }

  } else {

    ExecEnv::log().error("MutateGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Offset: {}, Size: {}",
                         genome_variant_ptr->genomeId(), contig, offset, region_size);
    return false;

  }

  return true;

}

