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
                                   std::shared_ptr<const GenomeDatabase> genome_db_ptr) {



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

  std::vector<std::shared_ptr<const AminoSequence>> amino_vector;
  for (auto summary : gene_summary_map) {

    for (auto mutant : summary.second.sequence_mutant_vec) {

      amino_vector.push_back(mutant);

    }

  }

  std::string amino_compare = AminoSequence::multipleCompare(amino_vector);
  ExecEnv::log().info("Comparison of all Amino sequences for Contig: {}, Gene: {}, Sequence: {}; \n{}",
                      contig, gene, sequence, amino_compare);



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

