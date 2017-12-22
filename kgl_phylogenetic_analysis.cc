//
// Created by kellerberrin on 17/11/17.
//

#include "kgl_phylogenetic_analysis.h"
#include "kgl_sequence_virtual_compare.h"

namespace kgl = kellerberrin::genome;




bool kgl::ApplicationAnalysis::writeMutantProteins(const std::string& fasta_file,
                                                   const std::string& sequence_name,
                                                   const ContigId_t& contig_id,
                                                   const FeatureIdent_t& gene_id,
                                                   const FeatureIdent_t& sequence_id,
                                                   const std::shared_ptr<const GenomeDatabase>& genome_db,
                                                   const std::shared_ptr<const GenomeVariant>& genome_variant) {


  std::vector<std::shared_ptr<AminoSequence>> mutant_sequence_vector;
  std::shared_ptr<AminoSequence> reference_sequence;

  if (genome_variant->mutantProteins(contig_id,
                                     gene_id,
                                     sequence_id,
                                     genome_db,
                                     reference_sequence,
                                     mutant_sequence_vector)) {

    std::vector<WriteFastaSequence> fasta_sequence_vec;
    size_t sequence_counter = 0;
    for (auto amino_sequence : mutant_sequence_vector) {

      ++sequence_counter;
      std::stringstream ss;
      ss << sequence_name << "_" << sequence_counter;
      WriteFastaSequence fasta_sequence;
      fasta_sequence.first = ss.str();
      fasta_sequence.second = amino_sequence;
      fasta_sequence_vec.push_back(fasta_sequence);

    }

    return ParseGffFasta().writeFastaFile(fasta_file, fasta_sequence_vec);

  } else {

    ExecEnv::log().warn("No valid sequence for contig: {}, gene: {}, sequence id: {} Fasta: {} not written",
                        contig_id, gene_id, sequence_id, fasta_file);
    return false;

  }

}


bool kgl::ApplicationAnalysis::readMutantProteins(const std::string& fasta_file,
                                                  const std::string& sequence_name,
                                                  const ContigId_t& contig_id,
                                                  const FeatureIdent_t& gene_id,
                                                  const FeatureIdent_t& sequence_id,
                                                  const std::shared_ptr<const GenomeDatabase>& genome_db,
                                                  const std::shared_ptr<const GenomeVariant>& genome_variant,
                                                  std::vector<std::string>& comparison_string_vector) {


  std::vector<std::shared_ptr<AminoSequence>> mutant_sequence_vector;
  std::shared_ptr<AminoSequence> reference_sequence;

  comparison_string_vector.clear();

  if (genome_variant->mutantProteins(contig_id,
                                     gene_id,
                                     sequence_id,
                                     genome_db,
                                     reference_sequence,
                                     mutant_sequence_vector)) {

    std::vector<ReadFastaSequence> fasta_sequence_vec;
    if (ParseGffFasta().readFastaFile(fasta_file, fasta_sequence_vec)) {

      for (auto sequence: fasta_sequence_vec) {

        if (sequence.first == sequence_name) {

          for (auto amino_sequence : mutant_sequence_vector) {

            StringAminoAcid fasta_amino_string(*sequence.second);
            std::shared_ptr<AminoSequence> fasta_amino_sequence(std::make_shared<AminoSequence>(fasta_amino_string));
            std::string comparison_string = fasta_amino_sequence->compareAminoSequences(*amino_sequence);
            comparison_string_vector.push_back(comparison_string);

          }

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


bool kgl::ApplicationAnalysis::compareMutantProteins(const ContigId_t& contig_id,
                                                     const FeatureIdent_t& gene_id,
                                                     const FeatureIdent_t& sequence_id,
                                                     const std::shared_ptr<const GenomeDatabase>& genome_db,
                                                     const std::shared_ptr<const GenomeVariant>& genome_variant,
                                                     std::vector<std::string>& comparison_string_vector) {

  std::vector<std::shared_ptr<AminoSequence>> mutant_sequence_vector;
  std::shared_ptr<AminoSequence> reference_sequence;

  comparison_string_vector.clear();

  if (genome_variant->mutantProteins(contig_id,
                                     gene_id,
                                     sequence_id,
                                     genome_db,
                                     reference_sequence,
                                     mutant_sequence_vector)) {

    for (auto amino_sequence : mutant_sequence_vector) {

      std::string comparison_string = reference_sequence->compareAminoSequences(*amino_sequence);
      comparison_string_vector.push_back(comparison_string);

    }

  } else {

    ExecEnv::log().warn("No valid sequence for contig: {}, gene: {}, sequence id: {}", contig_id, gene_id, sequence_id);
    return false;

  }

  return true;

}


bool kgl::PhylogeneticAnalysis::UPGMA(const std::string& newick_file,
                                      std::shared_ptr<const PopulationStatistics> population_stats_ptr) {

  UPGMAMatrix<const GenomeStatistics> upgma_matrix(population_stats_ptr->initUPGMA());

  upgma_matrix.calculateReduce();

  upgma_matrix.writeNewick(newick_file);

  return true;

}
