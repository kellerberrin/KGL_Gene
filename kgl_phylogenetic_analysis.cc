//
// Created by kellerberrin on 17/11/17.
//

#include "kgl_phylogenetic_analysis.h"
#include "kgl_sequence_virtual_compare.h"
#include "kgl_sequence_offset.h"

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

    std::string reference_name =  sequence_name + "_reference";
    WriteFastaSequence gff_sequence;
    gff_sequence.first = reference_name;
    gff_sequence.second = reference_sequence;
    fasta_sequence_vec.push_back(gff_sequence);

    for (auto amino_sequence : mutant_sequence_vector) {

      ++sequence_counter;
      std::stringstream ss;
      ss << sequence_name << "_mutant_" << sequence_counter;
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


bool kgl::ApplicationAnalysis::compareMutantCodingDNA(const ContigId_t& contig_id,
                                                     const FeatureIdent_t& gene_id,
                                                     const FeatureIdent_t& sequence_id,
                                                     const std::shared_ptr<const GenomeDatabase>& genome_db,
                                                     const std::shared_ptr<const GenomeVariant>& genome_variant,
                                                     std::vector<std::string>& comparison_string_vector) {

  std::vector<std::shared_ptr<DNA5SequenceCoding>> mutant_sequence_vector;
  std::shared_ptr<DNA5SequenceCoding> reference_sequence;

  comparison_string_vector.clear();

  if (genome_variant->mutantCodingDNA(contig_id,
                                     gene_id,
                                     sequence_id,
                                     genome_db,
                                     reference_sequence,
                                     mutant_sequence_vector)) {

    for (auto dna_sequence : mutant_sequence_vector) {

      std::string comparison_string = reference_sequence->compareDNA5Coding(*dna_sequence);
      comparison_string_vector.push_back(comparison_string);

    }

  } else {

    ExecEnv::log().warn("No valid sequence for contig: {}, gene: {}, sequence id: {}", contig_id, gene_id, sequence_id);
    return false;

  }

  return true;

}



bool kgl::ApplicationAnalysis::compareMutantRegions(const ContigId_t& contig_id,
                                                    ContigOffset_t region_offset,
                                                    ContigSize_t region_size,
                                                    StrandSense strand,
                                                    const std::shared_ptr<const GenomeDatabase>& genome_db,
                                                    const std::shared_ptr<const GenomeVariant>& genome_variant,
                                                    std::vector<std::string>& comparison_string_vector) {

  std::vector<std::shared_ptr<DNA5SequenceLinear>> mutant_sequence_vector;
  std::shared_ptr<DNA5SequenceLinear> reference_sequence;

  comparison_string_vector.clear();

  if (genome_variant->mutantRegion(contig_id,
                                   region_offset,
                                   region_size,
                                   genome_db,
                                   reference_sequence,
                                   mutant_sequence_vector)) {

    std::shared_ptr<DNA5SequenceCoding> reference_sequence_coding = SequenceOffset::codingSequence(reference_sequence, strand);
    for (auto dna_sequence : mutant_sequence_vector) {

      std::shared_ptr<DNA5SequenceCoding> mutant_sequence_coding = SequenceOffset::codingSequence(dna_sequence, strand);
      std::string comparison_string = reference_sequence_coding->compareDNA5Coding(*mutant_sequence_coding);
      comparison_string_vector.push_back(comparison_string);

    }

  } else {

    ExecEnv::log().warn("No valid sequence for contig: {}, offset: {}, size: {}", contig_id, region_offset, region_size);
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
                                             std::vector<std::string>& comparison_string_vector) {

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

  std::vector<std::shared_ptr<DNA5SequenceLinear>> mutant_sequence_vector;
  std::shared_ptr<DNA5SequenceLinear> reference_sequence;

  comparison_string_vector.clear();

  if (genome_variant->mutantRegion(contig_id,
                                   offset_5_prime,
                                   size_5_prime,
                                   genome_db,
                                   reference_sequence,
                                   mutant_sequence_vector)) {

    std::shared_ptr<DNA5SequenceCoding> reference_sequence_coding = SequenceOffset::codingSequence(reference_sequence, coding_sequence_ptr->strand());
    for (auto dna_sequence : mutant_sequence_vector) {

      std::shared_ptr<DNA5SequenceCoding> sequence_coding = SequenceOffset::codingSequence(dna_sequence, coding_sequence_ptr->strand());
      std::string comparison_string = reference_sequence_coding->compareDNA5Coding(*sequence_coding);
      comparison_string_vector.push_back(comparison_string);

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
                                             std::vector<std::string>& comparison_string_vector) {

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


  std::vector<std::shared_ptr<DNA5SequenceLinear>> mutant_sequence_vector;
  std::shared_ptr<DNA5SequenceLinear> reference_sequence;

  comparison_string_vector.clear();

  if (genome_variant->mutantRegion(contig_id,
                                   offset_3_prime,
                                   size_3_prime,
                                   genome_db,
                                   reference_sequence,
                                   mutant_sequence_vector)) {

    std::shared_ptr<DNA5SequenceCoding> reference_sequence_coding = SequenceOffset::codingSequence(reference_sequence, coding_sequence_ptr->strand());
    for (auto dna_sequence : mutant_sequence_vector) {

      std::shared_ptr<DNA5SequenceCoding> sequence_coding = SequenceOffset::codingSequence(dna_sequence, coding_sequence_ptr->strand());
      std::string comparison_string = reference_sequence_coding->compareDNA5Coding(*sequence_coding);
      comparison_string_vector.push_back(comparison_string);

    }

  } else {

    ExecEnv::log().warn("No valid 3 prime sequence for contig: {}, offset: {}, size: {}", contig_id, offset_3_prime, size_3_prime);
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
