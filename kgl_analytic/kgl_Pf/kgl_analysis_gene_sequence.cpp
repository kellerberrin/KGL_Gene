//
// Created by kellerberrin on 16/01/18.
//

#include "kgl_analysis_gene_sequence.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_sequence_complexity.h"
#include "kgl_seq_variant_db.h"
#include "kgl_seq_variant_filter.h"
#include "kgl_seq_coding.h"
#include "kgl_io_gff_fasta.h"

#include <memory>
#include <sstream>


namespace kgl = kellerberrin::genome;



bool kgl::GenomicSequence::translateContig(const GenomeId_t& genome_id,
                                           const ContigId_t& contig_id,
                                           const std::shared_ptr<const GenomeCollection>& genome_collection_ptr,
                                           const std::string& fasta_file_name) {
  // Get the genome object
  std::optional<std::shared_ptr<const GenomeReference>> genome_opt = genome_collection_ptr->getOptionalGenome(genome_id);
  if (not genome_opt) {

    ExecEnv::log().warn("Could not find Genome: {} in genome collection", genome_id);
    return false;
    
  }
  auto& genome_ref_ptr = genome_opt.value();

  auto contig_opt = genome_ref_ptr->getContigSequence(contig_id);
  if (not contig_opt) {

    ExecEnv::log().warn("Did not find Contig: {} in Genome: {}", contig_id, genome_id);
    return false;

  }
  auto& contig_ref_ptr = contig_opt.value();

  std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>> amino_fasta_vector;  
  for (const auto& [gene_id, gene_ptr] : contig_ref_ptr->getGeneMap()) {
    
    auto transcript_array_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);

    size_t coding_count = 0;
    for (const auto& [transcript_id, transcript_ptr] : transcript_array_ptr->getMap()) {

      ExecEnv::log().info("Sequence Id: {}", transcript_id);

      for (const auto& [cds_offset, cds_ptr] : transcript_ptr->getFeatureMap()) {

        ExecEnv::log().info("CDS Id: {}, Offset: {}", cds_ptr->id(), cds_offset);

        for (const auto& [attribute_key, attribute_value] : cds_ptr->getAttributes().getMap()) {

          ExecEnv::log().info("CDS Attribute Key: {} Value: {}", attribute_key, attribute_value);

        }

      }
    

      auto coding_dna_opt = contig_ref_ptr->codingSequence(transcript_ptr);
      if (coding_dna_opt) {

        DNA5SequenceCoding& coding_dna = coding_dna_opt.value();
        std::stringstream ss;
        ss << gene_ptr->id() << "_" << ++coding_count;
        AminoSequence peptide = transcript_ptr->contig()->getAminoSequence(coding_dna);
        std::pair<std::string, std::shared_ptr<AminoSequence>> fasta_entry(ss.str(), std::make_shared<AminoSequence>(std::move(peptide)));
        amino_fasta_vector.push_back(fasta_entry);

      } else {

        ExecEnv::log().error("Cannot generate a DNA coding sequence for Gene: {} in Genome: {}", gene_ptr->id(), genome_id);
        return false;

      } // Coding DNA 

    } // For all sequences

  }

  // Sequences are presented as a pair of a sequence name and an amino sequences.
  if (not GenomicMutation::writeMutantProteins(fasta_file_name, amino_fasta_vector)) {

    ExecEnv::log().warn("Problem writing peptide fasta for Contig: {} in file: {}", contig_id, fasta_file_name);
    return false;

  }

  return true;

}


bool kgl::GenomicSequence::mutateGene(const ContigId_t& contig,
                                      const FeatureIdent_t& gene,
                                      const FeatureIdent_t& sequence,
                                      const std::shared_ptr<const PopulationDB>& population_ptr,
                                      const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                      const std::string& fasta_filename) {



  GeneSummaryMap gene_summary_map;

  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    if (not mutateGenomeGene(contig, gene, sequence, genome_ptr, genome_ref_ptr, gene_summary_map)) {

      return false;

    }

  }

  // Get the contig_ref_ptr.
  std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_ref_ptr->getContigSequence(contig);
  if (not contig_opt) {

    ExecEnv::log().warn("Could not find contig_ref_ptr: {} in genome database", contig);
    return false;

  }
  auto& contig_ptr = contig_opt.value();

  for (auto summary : gene_summary_map) {

    auto sequence_validity = contig_ptr->checkValidProteinSequence(*summary.second.sequence_mutant);
    if (TranscriptionSequence::checkValidProtein(sequence_validity)) {

      ExecEnv::log().info("Multiple comparison Genome: {}, protein score: {}, 5Prime score: {}, 3Prime score: {}",
                          summary.first, summary.second.sequence_distance, summary.second.prime5_distance, summary.second.prime3_distance);

    } else {

      ExecEnv::log().info("Frame shift mutation Genome: {}, protein score: {}, 5Prime score: {}, 3Prime score: {}",
                          summary.first, summary.second.sequence_distance, summary.second.prime5_distance, summary.second.prime3_distance);


    }

  }

  // Write all the amino sequences as a fasta file.
  std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>> amino_fasta_vector;
  for (auto const& summary : gene_summary_map) {

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
  if (not GenomicMutation::writeMutantProteins(fasta_filename, amino_fasta_vector)) {

    ExecEnv::log().warn("Problem writing protein fasta file: {}", fasta_filename);

  }


  return true;

}



bool kgl::GenomicSequence::mutateGenomeGene(const ContigId_t& contig,
                                            const FeatureIdent_t& gene,
                                            const FeatureIdent_t& sequence,
                                            const std::shared_ptr<const GenomeDB>& genome_variant_ptr,
                                            const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                            GeneSummaryMap& gene_summary_map) {


  GeneSummary gene_summary;
  LevenshteinGlobal dna_distance_metric;
  LevenshteinGlobal amino_distance_metric;
  gene_summary.genome = genome_variant_ptr->genomeId();

  DNA5SequenceCoding prime5_reference;
  DNA5SequenceCoding prime5_mutant;

  if (GenomicMutation::compare5Prime(contig,
                                     gene,
                                     sequence,
                                     PRIME_REGION_SIZE,
                                     genome_ref_ptr,
                                     genome_variant_ptr,
                                     prime5_reference,
                                     prime5_mutant)) {



    gene_summary.prime5_distance = dna_distance_metric.coding_distance(prime5_reference, prime5_mutant);

    ExecEnv::log().info("5PRIME Genome: {}, Contig: {}, Gene: {}, Sequence: {} Levenshtein: {}, comparison:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.prime5_distance);

    gene_summary.prime5_reference = std::make_shared<DNA5SequenceCoding>(std::move(prime5_reference));
    gene_summary.prime5_mutant = std::make_shared<DNA5SequenceCoding>(std::move(prime5_mutant));


  } else {

    ExecEnv::log().error("Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                         genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }

  AminoSequence protein_reference;
  AminoSequence protein_mutant;

  if (GenomeMutation::mutantProteins(contig,
                                     gene,
                                     sequence,
                                     genome_ref_ptr,
                                     gene_summary.variant_map,
                                     protein_reference,
                                     protein_mutant)) {



    gene_summary.sequence_distance = amino_distance_metric.amino_distance( protein_reference, protein_mutant);
    ExecEnv::log().info("Genome: {}, Contig: {}, Gene: {}, Sequence: {} Levenshtein: {}, comparison:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.sequence_distance);

    gene_summary.sequence_ptr = std::make_shared<AminoSequence>(std::move(protein_reference));
    gene_summary.sequence_mutant = std::make_shared<AminoSequence>(std::move(protein_mutant));


  } else {

    ExecEnv::log().error("Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }

  DNA5SequenceCoding prime3_reference;
  DNA5SequenceCoding prime3_mutant;

  if (GenomicMutation::compare3Prime(contig,
                                     gene,
                                     sequence,
                                     PRIME_REGION_SIZE,
                                     genome_ref_ptr,
                                     genome_variant_ptr,
                                     prime3_reference,
                                     prime3_mutant)) {

    gene_summary.prime3_distance = dna_distance_metric.coding_distance(prime3_reference, prime3_mutant);
    ExecEnv::log().info("3PRIME Genome: {}, Contig: {}, Gene: {}, Sequence: {} Levenshtein: {}, comparison:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.prime3_distance);

    gene_summary.prime3_reference = std::make_shared<DNA5SequenceCoding>(std::move(prime3_reference));
    gene_summary.prime3_mutant = std::make_shared<DNA5SequenceCoding>(std::move(prime3_mutant));

  } else {

    ExecEnv::log().error("Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                         genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }

  auto result = gene_summary_map.insert(std::pair<GenomeId_t, GeneSummary>(gene_summary.genome, gene_summary));

  if (not result.second) {

    ExecEnv::log().error("Duplicate sequence; Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                         genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }


  return true;

}


bool kgl::GenomicSequence::mutateAllRegions(const std::string& file_name,
                                            ContigSize_t region_size,
                                            const std::shared_ptr<const LinearDNASequenceDistance>& dna_distance_metric,
                                            const std::shared_ptr<const PopulationDB>& pop_variant_ptr,
                                            const std::shared_ptr<const GenomeReference>& genome_db_ptr) {

  const char CSV_delimiter = ',';
  // open the file.
  std::fstream out_file(file_name, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  out_file << outputRegionHeader(CSV_delimiter) << '\n';

  for(auto const& genome_variant : pop_variant_ptr->getMap()) {

    ExecEnv::log().info("Processing genome: {}", genome_variant.first);
    ContigSize_t sequence_count = 0;

    for (auto contig : genome_db_ptr->getMap()) {

      for (ContigOffset_t offset = 0; offset < (contig.second->sequence().length() - region_size); offset += region_size) {

        out_file << outputGenomeRegion(CSV_delimiter, dna_distance_metric, contig.first, offset, region_size, genome_variant.second, genome_db_ptr) << '\n';
        ++sequence_count;

      }

    }

    ExecEnv::log().info("Genome: {} mutated: {} sequences.", genome_variant.first, sequence_count);

  }

  return out_file.good();

}


std::string kgl::GenomicSequence::outputRegionHeader(char delimiter) {

  std::stringstream ss;

  ss << "Genome" << delimiter;
  ss << "Contig" << delimiter;
  ss << "ContigSize" << delimiter;
  ss << "ContigOffset" << delimiter;
  ss << "RegionSize" << delimiter;
  ss << "Distance" << delimiter;
  ss << "CpG" << delimiter;
  std::vector<DNA5::Alphabet> symbols = DNA5::enumerateAlphabet();
  for (auto symbol : symbols) {

    ss << static_cast<char>(symbol) << delimiter;

  }

  return ss.str();

}


std::string kgl::GenomicSequence::outputGenomeRegion(char delimiter,
                                                     const std::shared_ptr<const LinearDNASequenceDistance>& dna_distance_metric,
                                                     const ContigId_t& contig_id,
                                                     const ContigOffset_t offset,
                                                     const ContigSize_t region_size,
                                                     const std::shared_ptr<const GenomeDB>& genome_db_ptr,
                                                     const std::shared_ptr<const GenomeReference>& genome_ref_ptr) {

  std::stringstream ss;

  auto contig_ref_opt = genome_ref_ptr->getContigSequence(contig_id);
  if (not contig_ref_opt) {

    ExecEnv::log().warn("Contig: {} not found for Genome: {}", contig_id, genome_ref_ptr->genomeId());
    return "<error>";

  }
  auto& contig_ref_ptr = contig_ref_opt.value();

  if (offset + region_size >= contig_ref_ptr->sequence().length()) {

    ExecEnv::log().error("Offset: {} + Region Size: {} exceed the Contig: {} size: {}",
                         offset, region_size, contig_id, contig_ref_ptr->sequence().length());
    return "<error>";

  }

  auto contig_db_opt = genome_db_ptr->getContig(contig_id);
  if (not contig_db_opt) {

    ExecEnv::log().warn("Contig: {} not found for Genome: {}", contig_id, genome_db_ptr->genomeId());
    return "<error>";

  }
  auto& contig_db_ptr = contig_db_opt.value();


  // Filter the mutating varinats.
  OpenRightUnsigned region_interval(offset, offset+region_size);
  SequenceVariantFilter seq_variant_filter(contig_db_ptr, region_interval);

  // And mutate the sequence.
  AdjustedSequence adjusted_sequence;
  if (not adjusted_sequence.updateSequence(contig_ref_ptr, seq_variant_filter)) {

    ExecEnv::log().warn("Problem mutating region DNA sequence for contig_id id: {}, interval: {}",
                        contig_db_ptr->contigId(), seq_variant_filter.sequenceInterval().toString());
    return "<error>";

  }

  auto dual_seq_opt = adjusted_sequence.moveSequenceClear();
  if (not dual_seq_opt) {

    ExecEnv::log().error("Unexpected error mutating Genome: {}, Contig: {}, Offset: {}, Size: {}",
                         genome_ref_ptr->genomeId(), contig_id, offset, region_size);
    return "<error>";

  }
  auto& [reference_sequence, mutant_sequence] = dual_seq_opt.value();

  double distance = static_cast<double>(dna_distance_metric->linear_distance(reference_sequence, mutant_sequence));

  ss << genome_db_ptr->genomeId() << delimiter;
  ss << contig_id << delimiter;
  ss << contig_ref_ptr->sequence().length() << delimiter;
  ss << offset << delimiter;
  ss << region_size << delimiter;
  ss << distance << delimiter;
  ss << SequenceComplexity::relativeCpGIslands(reference_sequence) << delimiter;
  std::vector<std::pair<DNA5::Alphabet, size_t>> symbol_count = reference_sequence.countSymbols();
  for (auto const& count : symbol_count) {

    ss << count.second << delimiter;

  }

  return ss.str();

}



bool kgl::GenomicSequence::mutateGenomeRegion(const GenomeId_t& genome,
                                              const ContigId_t& contig,
                                              ContigOffset_t offset,
                                              ContigSize_t region_size,
                                              const std::shared_ptr<const PopulationDB>& population_ptr,
                                              const std::shared_ptr<const GenomeReference>& genome_db_ptr,
                                              const std::string& fasta_file) {


  auto genome_opt = population_ptr->getGenome(genome);
  if (genome_opt) {

    return mutateGenomeRegion(contig, offset, region_size, genome_opt.value(), genome_db_ptr, fasta_file);

  } else {

    ExecEnv::log().error("Genome: {} not found", genome);
    return false;

  }

}


bool kgl::GenomicSequence::mutateGenomeRegion(const ContigId_t& contig_id,
                                              const ContigOffset_t offset,
                                              const ContigSize_t region_size,
                                              const std::shared_ptr<const GenomeDB>& genome_db_ptr,
                                              const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                              const std::string& fasta_file) {

  std::shared_ptr<const LinearDNASequenceDistance> dna_distance_metric(std::make_shared<const LevenshteinGlobal>());
  std::vector<WriteFastaSequence> dna_seq_vector;

  auto contig_db_opt = genome_db_ptr->getContig(contig_id);
  if (not contig_db_opt) {

    ExecEnv::log().warn("Contig: {} not found for Genome: {}", contig_id, genome_db_ptr->genomeId());
    return false;

  }
  const auto& contig_db_ptr = contig_db_opt.value();

  auto contig_ref_opt = genome_ref_ptr->getContigSequence(contig_id);
  if (not contig_db_opt) {

    ExecEnv::log().warn("Contig: {} not found for Genome: {}", contig_id, genome_ref_ptr->genomeId());
    return false;

  }
  const auto& contig_ref_ptr = contig_ref_opt.value();

  // Filter the mutating variants.
  OpenRightUnsigned region_interval(offset, offset+region_size);
  SequenceVariantFilter seq_variant_filter(contig_db_ptr, region_interval);

  // And mutate the sequence.
  AdjustedSequence adjusted_sequence;
  if (not adjusted_sequence.updateSequence(contig_ref_ptr, seq_variant_filter)) {

    ExecEnv::log().warn("Problem mutating region DNA sequence for contig_id id: {}, interval: {}",
                        contig_db_ptr->contigId(), seq_variant_filter.sequenceInterval().toString());
    return false;

  }

  auto dual_seq_opt = adjusted_sequence.moveSequenceClear();
  if (not dual_seq_opt) {

    ExecEnv::log().error("Unexpected error mutating Genome: {}, Contig: {}, Offset: {}, Size: {}",
                         genome_db_ptr->genomeId(), contig_id, offset, region_size);
    return false;

  }
  auto& [reference_sequence, mutant_sequence] = dual_seq_opt.value();

  DNA5SequenceCoding ref_reverse_complement = reference_sequence.codingSequence(StrandSense::REVERSE);

  CompareDistance_t score = dna_distance_metric->linear_distance(reference_sequence, mutant_sequence);

  ExecEnv::log().info("Genome: {}, Contig: {}, Offset: {} Size: {} score: {}",
                      genome_db_ptr->genomeId(), contig_id, offset, region_size, score);

  auto reference_sequence_ptr = std::make_shared<DNA5SequenceLinear>(std::move(reference_sequence));
  auto mutant_sequence_ptr = std::make_shared<DNA5SequenceLinear>(std::move(mutant_sequence));
  auto ref_reverse_complement_ptr = std::make_shared<DNA5SequenceCoding>(std::move(ref_reverse_complement));

  std::stringstream ref_fasta_ss;
  ref_fasta_ss << "reference_" << contig_id << "_+_" << offset << "_" << region_size;
  WriteFastaSequence fasta_reference(ref_fasta_ss.str(), "", reference_sequence_ptr);
  dna_seq_vector.push_back(fasta_reference);

  ref_fasta_ss.str("");
  ref_fasta_ss << "reference_" << contig_id << "_-_" << offset << "_" << region_size;
  WriteFastaSequence rev_fasta_reference(ref_fasta_ss.str(), "", ref_reverse_complement_ptr);
  dna_seq_vector.push_back(rev_fasta_reference);

  std::stringstream mutant_fasta_ss;
  mutant_fasta_ss << "mutant_" << contig_id << "_" << offset << "_" << region_size;

  WriteFastaSequence fasta_mutant(mutant_fasta_ss.str(), "", mutant_sequence_ptr);
  dna_seq_vector.push_back(fasta_mutant);

  if (not ParseFasta::writeFastaFile(fasta_file, dna_seq_vector)) {

    ExecEnv::log().error("Problem writing to fasta file: {}", fasta_file);

  }

  return true;

}

