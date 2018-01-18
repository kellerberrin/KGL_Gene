//
// Created by kellerberrin on 16/01/18.
//

#include "kgl_phylogenetic_gene.h"


namespace kgl = kellerberrin::genome;


bool kgl::GeneAnalysis::mutateGene(const ContigId_t& contig,
                                   const FeatureIdent_t& gene,
                                   const FeatureIdent_t& sequence,
                                   std::shared_ptr<const PopulationVariant> population_ptr,
                                   std::shared_ptr<const GenomeDatabase> genome_db_ptr) {


  for (auto genome : population_ptr->getMap()) {

    if (not mutateGenomeGene(contig, gene, sequence, genome.second, genome_db_ptr)) {

      return false;

    }

  }

  return true;

}


bool kgl::GeneAnalysis::mutateGenomeGene(const GenomeId_t& genome,
                                         const ContigId_t& contig,
                                         const FeatureIdent_t& gene,
                                         const FeatureIdent_t& sequence,
                                         std::shared_ptr<const PopulationVariant> population_ptr,
                                         std::shared_ptr<const GenomeDatabase> genome_db_ptr) {


  std::shared_ptr<const GenomeVariant> genome_variant_ptr;
  if (population_ptr->getGenomeVariant(genome, genome_variant_ptr)) {

    return mutateGenomeGene(contig, gene, sequence, genome_variant_ptr, genome_db_ptr);

  } else {

    ExecEnv::log().error("mutateGenomeGene(), Genome: {} not found", genome);
    return false;

  }

}


bool kgl::GeneAnalysis::mutateGenomeGene(const ContigId_t& contig,
                                         const FeatureIdent_t& gene,
                                         const FeatureIdent_t& sequence,
                                         std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                                         std::shared_ptr<const GenomeDatabase> genome_db_ptr) {




  std::vector<std::shared_ptr<AminoSequence>> mutant_sequence_vector;
  std::shared_ptr<AminoSequence> reference_sequence;
  if (genome_variant_ptr->mutantProteins(contig,
                                         gene,
                                         sequence,
                                         genome_db_ptr,
                                         reference_sequence,
                                         mutant_sequence_vector)) {


    for (auto mutant : mutant_sequence_vector) {

      CompareScore_t score;
      std::string comparison = reference_sequence->compareAminoSequences(mutant, score);
      ExecEnv::log().info("Genome: {}, Contig: {}, Gene: {}, Sequence: {} score: {}, comparison:\n{}",
                          genome_variant_ptr->genomeId(), contig, gene, sequence, score, comparison);

    }

  } else {

    ExecEnv::log().error("MutateGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
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

