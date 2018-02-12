//
// Created by kellerberrin on 17/11/17.
//

#include "kgl_upgma.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_sequence_offset.h"
#include "kgl_phylogenetic_gene.h"

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

bool kgl::ApplicationAnalysis::writeMutantDNA(const std::string& fasta_file,
                                              const std::vector<std::pair<std::string, std::shared_ptr<DNA5SequenceLinear>>>& dna_seq_vector) {

  std::vector<WriteFastaSequence> fasta_sequence_vec;

  for (auto amino_sequence : dna_seq_vector) {

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
  OffsetVariantMap variant_map;

  if (genome_variant->mutantRegion(contig_id,
                                   offset_5_prime,
                                   size_5_prime,
                                   genome_db,
                                   variant_map,
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
  OffsetVariantMap variant_map;

  if (genome_variant->mutantRegion(contig_id,
                                   offset_3_prime,
                                   size_3_prime,
                                   genome_db,
                                   variant_map,
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

  out_file << GeneAnalysis::outputRegionHeader(CSV_delimiter) << CSV_delimiter;
  out_file << outputSequenceHeader(CSV_delimiter) << '\n';

  for( auto genome_variant : pop_variant_ptr->getMap()) {

    ExecEnv::log().info("outputSequenceCSV(), Processing genome: {}", genome_variant.first);
    size_t sequence_count = 0;

    for (auto contig : genome_db->getMap()) {

      std::shared_ptr<const CodingSequence> previous_seq_ptr = nullptr;

      for (auto gene : contig.second->getGeneMap()) {

        const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene.second);
        for (auto sequence : coding_seq_ptr->getMap()) {

          ContigSize_t front_porch_size;
          ContigOffset_t front_porch_offset;
          if (not previous_seq_ptr) {

            front_porch_offset = 0;
            front_porch_size = sequence.second->start();

          } else {

            front_porch_offset = previous_seq_ptr->end();
            if (front_porch_offset <= sequence.second->start()) {

              front_porch_size = sequence.second->start() - front_porch_offset;

            } else {

              ExecEnv::log().info("Sequence: {} end offset: {} overlaps sequence: {} begin offset: {}",
                                  previous_seq_ptr->getCDSParent()->id(),
                                  previous_seq_ptr->end(),
                                  sequence.second->getCDSParent()->id(),
                                  sequence.second->start());
              front_porch_size = 0;

            }

          }

          out_file << GeneAnalysis::outputGenomeRegion(CSV_delimiter, contig.first, front_porch_offset, front_porch_size, genome_variant.second, genome_db);
          out_file << CSV_delimiter;
          out_file << outputSequence(CSV_delimiter, sequence.second, genome_db, genome_variant.second);
          out_file << '\n';
          ++sequence_count;
          previous_seq_ptr = sequence.second;

        }

      }

    }

    ExecEnv::log().info("outputSequenceCSV(), Genome: {} mutated: {} sequences.", genome_variant.first, sequence_count);

  }

  return out_file.good();

}


std::string kgl::ApplicationAnalysis::outputSequenceHeader(char delimiter) {

  std::stringstream ss;

  ss << "Genome" << delimiter;
  ss << "Contig" << delimiter;
  ss << "ContigLength" << delimiter;
  ss << "Sequence" << delimiter;
  ss << "ContigOffset" << delimiter;
  ss << "DNASize" << delimiter;
  ss << "SequenceGC" << delimiter;
  ss << "Error" << delimiter;
  ss << "ValidReference" << delimiter;
  ss << "AllPaths" << delimiter;
  ss << "ValidPaths" << delimiter;
  ss << "Score" << delimiter;
  ss << "Strand" << delimiter;
  ss << "DNAScore" << delimiter;
  ss << "Symbolic" << delimiter;
  ss << "AltSymbolic" << delimiter;
  ss << "Description";

  return ss.str();

}


std::string kgl::ApplicationAnalysis::outputSequence(char delimiter,
                                                      std::shared_ptr<const CodingSequence> coding_sequence,
                                                      std::shared_ptr<const GenomeDatabase> genome_db,
                                                      std::shared_ptr<const GenomeVariant> genome_variant) {

  std::string genome_id = genome_variant->genomeId();
  std::shared_ptr<const ContigFeatures> contig_ptr = coding_sequence->getGene()->contig();
  std::string gene_id = coding_sequence->getGene()->id();
  std::string sequence_id = coding_sequence->getCDSParent()->id();
  ContigOffset_t sequence_offset = coding_sequence->start();
  char strand = static_cast<char>(coding_sequence->strand());

  bool error_flag = true;
  size_t mutant_paths = 0;
  size_t valid_paths = 0;
  double average_score = 0;
  double proportion_GC = 0;
  double average_DNA_score = 0;
  bool valid_reference = false;
  OffsetVariantMap variant_map;

  std::shared_ptr<DNA5SequenceCoding> reference_sequence;
  std::vector<std::shared_ptr<DNA5SequenceCoding>> mutant_sequence_vector;
  if (genome_variant->mutantCodingDNA( contig_ptr->contigId(),
                                       gene_id,
                                       sequence_id,
                                       genome_db,
                                       variant_map,
                                       reference_sequence,
                                       mutant_sequence_vector)) {

    for (auto mutant : mutant_sequence_vector) {

      CompareScore_t DNA_score;
      DNA_score = reference_sequence->compareMyerHirschberg(mutant);
      average_DNA_score += static_cast<double>(DNA_score);

    }

    if (reference_sequence->length() > 0) {

      proportion_GC = static_cast<double>(reference_sequence->countGC()) / static_cast<double>(reference_sequence->length());

    }

  }



  if (mutant_sequence_vector.size()) {

    average_DNA_score = average_DNA_score/ static_cast<double>(mutant_sequence_vector.size());

  }

  std::shared_ptr<AminoSequence> amino_reference_seq;
  std::vector<std::shared_ptr<AminoSequence>> amino_mutant_vec;
  if (genome_variant->mutantProteins(contig_ptr->contigId(),
                                     gene_id,
                                     sequence_id,
                                     genome_db,
                                     variant_map,
                                     amino_reference_seq,
                                     amino_mutant_vec)) {

    error_flag = false;
    valid_reference = contig_ptr->verifyProteinSequence(amino_reference_seq);
    for (auto mutant : amino_mutant_vec) {

      CompareScore_t amino_score;
      amino_reference_seq->compareAminoSequences(mutant, amino_score);
      if (contig_ptr->verifyProteinSequence(mutant)) {

        average_score += static_cast<double>(amino_score);
        ++valid_paths;

      }
      ++mutant_paths;

    }

    if (valid_paths > 0) {

      average_score = average_score / static_cast<double>(valid_paths);

    } else {

      average_score = 0;

    }


  } else {

    ExecEnv::log().error("Problem mutating contig: {}, sequence: {}", contig_ptr->contigId(), sequence_id);

  }


  std::stringstream ss;

  ss << genome_id << delimiter;
  ss << contig_ptr->contigId() << delimiter;
  ss << contig_ptr->sequence().length() << delimiter;
  ss << sequence_id << delimiter;
  ss << sequence_offset << delimiter;
  ss << coding_sequence->codingNucleotides() << delimiter;
  ss << proportion_GC << delimiter;
  ss << error_flag << delimiter;
  ss << valid_reference << delimiter;
  ss << mutant_paths << delimiter;
  ss << valid_paths << delimiter;
  ss << average_score << delimiter;
  ss << strand << delimiter;
  ss << average_DNA_score << delimiter;

  std::shared_ptr<const OntologyRecord> gene_ontology_ptr;
  if (genome_db->geneOntology().getGafFeatureVector(gene_id, gene_ontology_ptr)) {

    ss << gene_ontology_ptr->symbolicReference() << delimiter;
    ss << gene_ontology_ptr->altSymbolicReference() << delimiter;
    ss << "\"" << gene_ontology_ptr->description() << "\"";

  } else {

    ss << "<NULL>" << delimiter;
    ss << "<NULL>" << delimiter;
    std::vector<std::string> description_vec;
    coding_sequence->getGene()->getAttributes().getDescription(description_vec);
    for (const auto &description : description_vec) {

      ss << "\"" << description << "\"";

    }

  }


  return ss.str();

}
