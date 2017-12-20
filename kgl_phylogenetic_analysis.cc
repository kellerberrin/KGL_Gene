//
// Created by kellerberrin on 17/11/17.
//

#include "kgl_phylogenetic_analysis.h"

namespace kgl = kellerberrin::genome;




bool kgl::ApplicationAnalysis::writeMutantProtein(const std::string& fasta_file,
                                                  const std::string& sequence_name,
                                                  const ContigId_t& contig_id,
                                                  const FeatureIdent_t& gene_id,
                                                  const FeatureIdent_t& sequence_id,
                                                  const std::shared_ptr<const GenomeDatabase>& genome_db,
                                                  const std::shared_ptr<const GenomeVariant>& genome_variant) {


  std::shared_ptr<AminoSequence> amino_sequence;
  if (genome_variant->mutantProtein(contig_id, gene_id, sequence_id, genome_db, amino_sequence)) {

    WriteFastaSequence fasta_sequence;
    fasta_sequence.first = sequence_name;
    fasta_sequence.second = amino_sequence;
    std::vector<WriteFastaSequence> fasta_sequence_vec;
    fasta_sequence_vec.push_back(fasta_sequence);
    return ParseGffFasta().writeFastaFile(fasta_file, fasta_sequence_vec);

  } else {

    ExecEnv::log().warn("No valid sequence for contig: {}, gene: {}, sequence id: {} Fasta: {} not written",
                        contig_id, gene_id, sequence_id, fasta_file);
    return false;

  }

}


bool kgl::ApplicationAnalysis::readMutantProtein(const std::string& fasta_file,
                                                 const std::string& sequence_name,
                                                 const ContigId_t& contig_id,
                                                 const FeatureIdent_t& gene_id,
                                                 const FeatureIdent_t& sequence_id,
                                                 const std::shared_ptr<const GenomeDatabase>& genome_db,
                                                 const std::shared_ptr<const GenomeVariant>& genome_variant,
                                                 std::string& comparison_string) {


  std::shared_ptr<AminoSequence> amino_sequence;
  if (genome_variant->mutantProtein(contig_id, gene_id, sequence_id, genome_db, amino_sequence)) {

    std::vector<ReadFastaSequence> fasta_sequence_vec;
    if (ParseGffFasta().readFastaFile(fasta_file, fasta_sequence_vec)) {

      for (auto sequence: fasta_sequence_vec) {

        if (sequence.first == sequence_name) {

          StringAminoAcid fasta_amino_string(*sequence.second);
          std::shared_ptr<AminoSequence> fasta_amino_sequence(std::make_shared<AminoSequence>(fasta_amino_string));
          comparison_string = fasta_amino_sequence->compareSequences(amino_sequence);

          return true;

        }

      }

      ExecEnv::log().warn("No sequence: {} found in Fasta file: {}", sequence_name, fasta_file);
      return false;


    } else {

      ExecEnv::log().warn("Unable read Fasta file: {}", fasta_file);
      return false;

    }

  } else {

    ExecEnv::log().warn("No valid sequence for contig: {}, gene: {}, sequence id: {}",
                        contig_id, gene_id, sequence_id, fasta_file);
    return false;

  }

}



bool kgl::PhylogeneticAnalysis::UPGMA(const std::string& newick_file,
                                      std::shared_ptr<const PopulationStatistics> population_stats_ptr) {

  UPGMAMatrix<const GenomeStatistics> upgma_matrix(population_stats_ptr->initUPGMA());

  upgma_matrix.calculateReduce();

  upgma_matrix.writeNewick(newick_file);

  return true;

}
