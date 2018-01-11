//
// Created by kellerberrin on 17/11/17.
//

#include "kgl_phylogenetic_analysis.h"
#include "kgl_sequence_virtual_compare.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


bool kgl::ApplicationAnalysis::writeMutantProteins(const std::string& fasta_file,
                                                   const std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>>& amino_seq_vector) {

  std::vector<WriteFastaSequence> fasta_sequence_vec;

  for (auto amino_sequence : amino_seq_vector) {

    WriteFastaSequence fasta_sequence;
    fasta_sequence.first = amino_sequence.first;
    fasta_sequence.second = amino_sequence.second;
    fasta_sequence_vec.push_back(fasta_sequence);

  }

  return ParseGffFasta().writeFastaFile(fasta_file, fasta_sequence_vec);


}


bool kgl::ApplicationAnalysis::readFastaProteins(const std::string& fasta_file,
                                                 std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>>& amino_seq_vector) {


  std::vector<std::shared_ptr<AminoSequence>> mutant_sequence_vector;
  std::shared_ptr<AminoSequence> reference_sequence;

  std::vector<ReadFastaSequence> fasta_sequence_vec;
  if (ParseGffFasta().readFastaFile(fasta_file, fasta_sequence_vec)) {

    for (auto sequence: fasta_sequence_vec) {

      StringAminoAcid fasta_amino_string(*sequence.second);
      std::shared_ptr<AminoSequence> fasta_amino_seq_ptr(std::make_shared<AminoSequence>(fasta_amino_string));
      std::pair<std::string, std::shared_ptr<AminoSequence>> insert_pair(sequence.first, fasta_amino_seq_ptr);
      amino_seq_vector.push_back(insert_pair);

    }

  } else {

    ExecEnv::log().warn("Unable read Fasta file: {}", fasta_file);
    return false;

  }

  return true;

}




bool kgl::ApplicationAnalysis::compare5Prime(const ContigId_t& contig_id,
                                             const FeatureIdent_t& gene_id,
                                             const FeatureIdent_t& sequence_id,
                                             ContigSize_t region_size,
                                             const std::shared_ptr<const GenomeDatabase>& genome_db,
                                             const std::shared_ptr<const GenomeVariant>& genome_variant,
                                             std::shared_ptr<DNA5SequenceCoding>& reference_sequence,
                                             std::vector<std::shared_ptr<DNA5SequenceCoding>>& mutant_sequence_vector) {

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("compare5Prime(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the coding sequence.
  std::shared_ptr<const CodingSequence> coding_sequence_ptr;
  if (not contig_ptr->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("compare5Prime(), Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }

  ContigOffset_t offset_5_prime;
  ContigSize_t size_5_prime;
  coding_sequence_ptr->prime_5_region(region_size, offset_5_prime, size_5_prime);

  ExecEnv::log().info("Analyzing the 5 Prime Region for Sequence {} at offset {}, size {}, 5 prime {}",
                      coding_sequence_ptr->getCDSParent()->id(), offset_5_prime, size_5_prime, coding_sequence_ptr->prime_5());

  std::vector<std::shared_ptr<DNA5SequenceLinear>> linear_mutant_sequence_vector;
  std::shared_ptr<DNA5SequenceLinear> linear_reference_sequence;

  if (genome_variant->mutantRegion(contig_id,
                                   offset_5_prime,
                                   size_5_prime,
                                   genome_db,
                                   linear_reference_sequence,
                                   linear_mutant_sequence_vector)) {

    reference_sequence = SequenceOffset::codingSequence(linear_reference_sequence, coding_sequence_ptr->strand());
    for (auto dna_sequence : linear_mutant_sequence_vector) {

      std::shared_ptr<DNA5SequenceCoding> sequence_coding = SequenceOffset::codingSequence(dna_sequence, coding_sequence_ptr->strand());
      mutant_sequence_vector.push_back(sequence_coding);

    }

  } else {

    ExecEnv::log().warn("No valid 5 prime sequence for contig: {}, offset: {}, size: {}", contig_id, offset_5_prime, size_5_prime);
    return false;

  }

  return true;

}


bool kgl::ApplicationAnalysis::compare3Prime(const ContigId_t& contig_id,
                                             const FeatureIdent_t& gene_id,
                                             const FeatureIdent_t& sequence_id,
                                             ContigSize_t region_size,
                                             const std::shared_ptr<const GenomeDatabase>& genome_db,
                                             const std::shared_ptr<const GenomeVariant>& genome_variant,
                                             std::shared_ptr<DNA5SequenceCoding>& reference_sequence,
                                             std::vector<std::shared_ptr<DNA5SequenceCoding>>& mutant_sequence_vector) {

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("compare5Prime(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the coding sequence.
  std::shared_ptr<const CodingSequence> coding_sequence_ptr;
  if (not contig_ptr->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("compare5Prime(), Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }

  ContigOffset_t offset_3_prime;
  ContigSize_t size_3_prime;
  coding_sequence_ptr->prime_3_region(region_size, offset_3_prime, size_3_prime);

  ExecEnv::log().info("Analyzing the 3 Prime Region for Sequence {} offset at {}, size {}, 3 prime offset {}",
                      coding_sequence_ptr->getCDSParent()->id(), offset_3_prime, size_3_prime, coding_sequence_ptr->prime_3());


  std::vector<std::shared_ptr<DNA5SequenceLinear>> linear_mutant_sequence_vector;
  std::shared_ptr<DNA5SequenceLinear> linear_reference_sequence;

  if (genome_variant->mutantRegion(contig_id,
                                   offset_3_prime,
                                   size_3_prime,
                                   genome_db,
                                   linear_reference_sequence,
                                   linear_mutant_sequence_vector)) {

    reference_sequence = SequenceOffset::codingSequence(linear_reference_sequence, coding_sequence_ptr->strand());
    for (auto dna_sequence : linear_mutant_sequence_vector) {

      std::shared_ptr<DNA5SequenceCoding> sequence_coding = SequenceOffset::codingSequence(dna_sequence, coding_sequence_ptr->strand());
      mutant_sequence_vector.push_back(sequence_coding);

    }

  } else {

    ExecEnv::log().warn("No valid 3 prime sequence for contig: {}, offset: {}, size: {}", contig_id, offset_3_prime, size_3_prime);
    return false;

  }

  return true;

}


bool kgl::ApplicationAnalysis::outputSequenceCSV(const std::string &file_name,
                                                 std::shared_ptr<const GenomeDatabase> genome_db,
                                                 std::shared_ptr<const PopulationVariant> pop_variant_ptr) {

  const char CSV_delimiter = ',';
  // open the file.
  std::fstream out_file(file_name, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  out_file << outputSequenceHeader(CSV_delimiter);

  for( auto genome_variant : pop_variant_ptr->getMap()) {

    for (auto contig : genome_db->getMap()) {

      for (auto gene : contig.second->getGeneMap()) {

        const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene.second);
        for (auto sequence : coding_seq_ptr->getMap()) {

          out_file << outputSequence(CSV_delimiter, sequence.second, genome_db, genome_variant.second);

        }

      }

    }

  }

  return out_file.good();

}


std::string kgl::ApplicationAnalysis::outputSequenceHeader(char delimiter) {

  std::stringstream ss;

  ss << "Genome" << delimiter;
  ss << "Contig" << delimiter;
  ss << "Sequence" << delimiter;
  ss << "Size(DNA)" << delimiter;
  ss << "Error" << delimiter;
  ss << "Paths" << delimiter;
  ss << "FrameShift" << delimiter;
  ss << "Score" << delimiter;
  ss << "Description" << '\n';

  return ss.str();

}


std::string kgl::ApplicationAnalysis::outputSequence(char delimiter,
                                                      std::shared_ptr<const CodingSequence> coding_sequence,
                                                      std::shared_ptr<const GenomeDatabase> genome_db,
                                                      std::shared_ptr<const GenomeVariant> genome_variant) {

  std::string genome_id = genome_variant->genomeId();
  std::string contig_id = coding_sequence->getGene()->contig()->contigId();
  std::string gene_id = coding_sequence->getGene()->id();
  std::string sequence_id = coding_sequence->getCDSParent()->id();
  std::vector<std::pair<std::string,std::string>> description_vec;
  coding_sequence->getGene()->getAttributes().getAllAttributes(description_vec);

  std::shared_ptr<AminoSequence> amino_reference_seq;
  std::vector<std::shared_ptr<AminoSequence>> amino_mutant_vec;
  bool frame_shift_flag;
  bool error_flag = true;
  size_t mutant_paths = 0;
  double average_score = 0;
  if (genome_variant->mutantProteins(contig_id,
                                     gene_id,
                                     sequence_id,
                                     genome_db,
                                     frame_shift_flag,
                                     amino_reference_seq,
                                     amino_mutant_vec)) {

    error_flag = false;
    for (auto mutant : amino_mutant_vec) {

      CompareScore_t amino_score;
      amino_reference_seq->compareAminoSequences(mutant, amino_score);
      average_score += static_cast<double>(amino_score);
      ++mutant_paths;

    }

    if (mutant_paths > 0) {

      average_score = average_score / static_cast<double>(mutant_paths);

    } else {

      average_score = 0;

    }


  } else {

    ExecEnv::log().error("Problem mutating contig: {}, sequence: {}", contig_id, sequence_id);

  }


  std::stringstream ss;

  ss << genome_id << delimiter;
  ss << contig_id << delimiter;
  ss << sequence_id << delimiter;
  ss << coding_sequence->codingNucleotides() << delimiter;
  ss << error_flag << delimiter;
  ss << mutant_paths << delimiter;
  ss << frame_shift_flag << delimiter;
  ss << average_score << delimiter;

  for (const auto& description : description_vec) {

    if (description.first == Attributes::DESCRIPTION_KEY) {

      ss << description.second;

    }

  }

  ss << '\n';

  return ss.str();

}


bool kgl::PhylogeneticAnalysis::UPGMA(const std::string& newick_file,
                                      std::shared_ptr<const PopulationStatistics> population_stats_ptr) {

  UPGMAMatrix<const GenomeStatistics> upgma_matrix(population_stats_ptr->initUPGMA());

  upgma_matrix.calculateReduce();

  upgma_matrix.writeNewick(newick_file);

  return true;

}
