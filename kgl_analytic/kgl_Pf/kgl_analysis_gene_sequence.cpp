//
// Created by kellerberrin on 16/01/18.
//

#include "kgl_sequence_offset.h"
#include "kgl_analysis_gene_sequence.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_sequence_complexity.h"
#include "kgl_variant_db_mutation.h"

#include <memory>
#include <sstream>


namespace kgl = kellerberrin::genome;



bool kgl::GenomicSequence::translateContig(const GenomeId_t& genome_id,
                                           const ContigId_t& contig_id,
                                           const std::shared_ptr<const GenomeCollection>& genomes,
                                           const std::string& fasta_file_name) {
  // Get the genome object
  std::optional<std::shared_ptr<const GenomeReference>> genome_opt = genomes->getOptionalGenome(genome_id);
  if (not genome_opt) {

    ExecEnv::log().warn("Could not find Genome: {} in genome collection", genome_id);
    return false;
    
  }
  std::shared_ptr<const GenomeReference> genome_db_ptr = genome_opt.value();

  // Look through the contigs for the contig.
  bool found = false;
  std::shared_ptr<const ContigReference> contig_ptr;
  for (const auto& contig : genome_db_ptr->getMap() ) {

    contig_ptr = contig.second;
    if (contig_ptr->contigId() == contig_id) {

      found = true;
      break;   

    }

  }

  if (not found or not contig_ptr) {

    ExecEnv::log().warn("Did not find Contig: {} in Genome: {}", contig_id, genome_id);
    return false;

  }

  std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>> amino_fasta_vector;  
  for (const auto& gene : contig_ptr->getGeneMap()) {
    
    std::shared_ptr<const kgl::CodingSequenceArray> coding_sequence_array_ptr;
    coding_sequence_array_ptr = GeneFeature::getCodingSequences(gene.second); 

    size_t coding_count = 0;
    for (const auto& coding_sequence : coding_sequence_array_ptr->getMap()) {

      ExecEnv::log().info("Sequence Id", coding_sequence.first);

      for (const auto& cds : coding_sequence.second->getSortedCDS()) {

      ExecEnv::log().info("CDS Id", cds.first);

        for (const auto& attribute : cds.second->getAttributes().getMap()) {

          ExecEnv::log().info("CDS Attribute Key: {} Value: {}", attribute.first, attribute.second);

        }

      }
    
      DNA5SequenceCoding coding_dna;
      if (contig_ptr->getDNA5SequenceCoding(coding_sequence.second, coding_dna)) {

        std::stringstream ss;
        ss << gene.second->id() << "_" << ++coding_count;
        AminoSequence peptide = coding_sequence.second->contig()->getAminoSequence(coding_dna);
        std::pair<std::string, std::shared_ptr<AminoSequence>> fasta_entry(ss.str(), std::make_shared<AminoSequence>(std::move(peptide)));
        amino_fasta_vector.push_back(fasta_entry);

      } else {

        ExecEnv::log().error("Cannot generate a DNA coding sequence for Gene: {} in Genome: {}", gene.second->id(), genome_id);
        return false;

      } // Coding DNA 

    } // For all sequences

  }

  // Sequences are presented as a pair of a sequence name and an amino sequences.
  if (not ApplicationAnalysis::writeMutantProteins(fasta_file_name, amino_fasta_vector)) {

    ExecEnv::log().warn("Problem writing peptide fasta for Contig: {} in file: {}", contig_id, fasta_file_name);
    return false;

  }

  return true;

}


bool kgl::GenomicSequence::translateGene(const GenomeId_t& genome_id,
                                         const FeatureIdent_t& gene_id,
                                         const std::shared_ptr<const GenomeCollection>& genomes,
                                         const std::string& fasta_file_name) {

  // Get the genome object
  std::optional<std::shared_ptr<const GenomeReference>> genome_opt = genomes->getOptionalGenome(genome_id);
  if (not genome_opt) {

    ExecEnv::log().warn("Could not find Genome: {} in genome collection", genome_id);
    return false;
    
  }
  std::shared_ptr<const GenomeReference> genome_db_ptr = genome_opt.value();

  // Look through the contigs for the gene.
  bool found = false;
  std::shared_ptr<const ContigReference> contig_ptr;
  std::vector<std::shared_ptr<const Feature>> feature_ptr_vec;
  for (const auto& contig : genome_db_ptr->getMap() ) {

    contig_ptr = contig.second;
    if (contig_ptr->findFeatureId(gene_id, feature_ptr_vec)) {

      found = true;
      break;   

    }

  }

  if (not found or feature_ptr_vec.empty() or not contig_ptr) {

    ExecEnv::log().warn("Did not find Gene: {} in Genome: {}", gene_id, genome_id);
    return false;

  }

  std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>> amino_fasta_vector;
  for (const auto& feature_ptr : feature_ptr_vec) {

    
    if (feature_ptr->isGene()) {

      std::shared_ptr<const kgl::CodingSequenceArray> coding_sequence_array_ptr;
      coding_sequence_array_ptr = GeneFeature::getCodingSequences(std::static_pointer_cast<const GeneFeature>(feature_ptr)); 

      size_t coding_count = 0;
      for (const auto& coding_sequence : coding_sequence_array_ptr->getMap()) {

        DNA5SequenceCoding coding_dna;
        if (contig_ptr->getDNA5SequenceCoding(coding_sequence.second, coding_dna)) {

          AminoSequence peptide_ptr = coding_sequence.second->contig()->getAminoSequence(coding_dna);
          std::stringstream ss;
          ss << gene_id << "_" << ++coding_count;
          std::pair<std::string, std::shared_ptr<AminoSequence>> fasta_entry(ss.str(), std::make_shared<AminoSequence>(std::move(peptide_ptr)));
          amino_fasta_vector.push_back(fasta_entry);

        } else {

          ExecEnv::log().error("Cannot generate a DNA coding sequence for Gene: {} in Genome: {}", gene_id, genome_id);
          return false;

        }

      }

    }

  }

  // Sequences are presented as a pair of a sequence name and an amino sequences.
  if (not ApplicationAnalysis::writeMutantProteins(fasta_file_name, amino_fasta_vector)) {

    ExecEnv::log().warn("Problem writing peptide fasta for Gene: {} in file: {}", gene_id, fasta_file_name);
    return false;

  }

  return true;

}



bool kgl::GenomicSequence::mutateGene(const ContigId_t& contig,
                                      const FeatureIdent_t& gene,
                                      const FeatureIdent_t& sequence,
                                      const std::shared_ptr<const PopulationDB>& population_ptr,
                                      const std::shared_ptr<const GenomeReference>& genome_db_ptr,
                                      const std::string& fasta_filename) {



  GeneSummaryMap gene_summary_map;

  for (auto const& genome : population_ptr->getMap()) {

    if (not mutateGenomeGene(contig, gene, sequence, genome.second, genome_db_ptr, gene_summary_map)) {

      return false;

    }

  }

  // Get the contig.
  std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db_ptr->getContigSequence(contig);
  if (not contig_opt) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig);
    return false;

  }

  for (auto summary : gene_summary_map) {

    if (contig_opt.value()->verifyProteinSequence(*summary.second.sequence_mutant)) {

      ExecEnv::log().info("Multiple comparison Genome: {}, protein score: {}, 5Prime score: {}, 3Prime score: {}",
                          summary.first, summary.second.sequence_distance, summary.second.prime5_distance, summary.second.prime3_distance);

    } else {

      ExecEnv::log().info("Frame shift mutation Genome: {}, protein score: {}, 5Prime score: {}, 3Prime score: {}",
                          summary.first, summary.second.sequence_distance, summary.second.prime5_distance, summary.second.prime3_distance);


    }

  }

  // Write all the amino sequences as a fasta file.
  std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>> amino_fasta_vector;
  for (auto summary : gene_summary_map) {

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
  if (not ApplicationAnalysis::writeMutantProteins(fasta_filename, amino_fasta_vector)) {

    ExecEnv::log().warn("mutantProtein(), Problem writing protein fasta file: {}", fasta_filename);

  }


  return true;

}



bool kgl::GenomicSequence::mutateGenomeGene(const ContigId_t& contig,
                                            const FeatureIdent_t& gene,
                                            const FeatureIdent_t& sequence,
                                            const std::shared_ptr<const GenomeDB>& genome_variant_ptr,
                                            const std::shared_ptr<const GenomeReference>& genome_db_ptr,
                                            GeneSummaryMap& gene_summary_map) {


  GeneSummary gene_summary;
  LevenshteinGlobal dna_distance_metric;
  LevenshteinGlobal amino_distance_metric;
  gene_summary.genome = genome_variant_ptr->genomeId();

  DNA5SequenceCoding prime5_reference;
  DNA5SequenceCoding prime5_mutant;

  if (ApplicationAnalysis::compare5Prime(contig,
                                         gene,
                                         sequence,
                                         PRIME_REGION_SIZE,
                                         genome_db_ptr,
                                         genome_variant_ptr,
                                         prime5_reference,
                                         prime5_mutant)) {



    gene_summary.prime5_distance = dna_distance_metric.coding_distance(prime5_reference, prime5_mutant);

    ExecEnv::log().info("5PRIME Genome: {}, Contig: {}, Gene: {}, Sequence: {} Levenshtein: {}, comparison:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.prime5_distance);

    gene_summary.prime5_reference = std::make_shared<DNA5SequenceCoding>(std::move(prime5_reference));
    gene_summary.prime5_mutant = std::make_shared<DNA5SequenceCoding>(std::move(prime5_mutant));


  } else {

    ExecEnv::log().error("MutateGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                         genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }

  AminoSequence protein_reference;
  AminoSequence protein_mutant;

  if (GenomeMutation::mutantProteins( contig,
                                      gene,
                                      sequence,
                                      genome_db_ptr,
                                      gene_summary.variant_map,
                                      protein_reference,
                                      protein_mutant)) {



    gene_summary.sequence_distance = amino_distance_metric.amino_distance( protein_reference, protein_mutant);
    ExecEnv::log().info("Genome: {}, Contig: {}, Gene: {}, Sequence: {} Levenshtein: {}, comparison:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.sequence_distance);

    gene_summary.sequence_ptr = std::make_shared<AminoSequence>(std::move(protein_reference));
    gene_summary.sequence_mutant = std::make_shared<AminoSequence>(std::move(protein_mutant));


  } else {

    ExecEnv::log().error("MutateGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Gene: {}, Sequence: {}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence);
    return false;

  }

  DNA5SequenceCoding prime3_reference;
  DNA5SequenceCoding prime3_mutant;

  if (ApplicationAnalysis::compare3Prime(contig,
                                         gene,
                                         sequence,
                                         PRIME_REGION_SIZE,
                                         genome_db_ptr,
                                         genome_variant_ptr,
                                         prime3_reference,
                                         prime3_mutant)) {

    gene_summary.prime3_distance = dna_distance_metric.coding_distance(prime3_reference, prime3_mutant);
    ExecEnv::log().info("3PRIME Genome: {}, Contig: {}, Gene: {}, Sequence: {} Levenshtein: {}, comparison:\n{}",
                        genome_variant_ptr->genomeId(), contig, gene, sequence, gene_summary.prime3_distance);

    gene_summary.prime3_reference = std::make_shared<DNA5SequenceCoding>(std::move(prime3_reference));
    gene_summary.prime3_mutant = std::make_shared<DNA5SequenceCoding>(std::move(prime3_mutant));

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

    ExecEnv::log().info("outputSequenceCSV(), Processing genome: {}", genome_variant.first);
    ContigSize_t sequence_count = 0;

    for (auto contig : genome_db_ptr->getMap()) {

      for (ContigOffset_t offset = 0; offset < (contig.second->sequence().length() - region_size); offset += region_size) {

        out_file << outputGenomeRegion(CSV_delimiter, dna_distance_metric, contig.first, offset, region_size, genome_variant.second, genome_db_ptr) << '\n';
        ++sequence_count;

      }

    }

    ExecEnv::log().info("outputSequenceCSV(), Genome: {} mutated: {} sequences.", genome_variant.first, sequence_count);

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
                                                     const std::shared_ptr<const GenomeDB>& genome_variant_ptr,
                                                     const std::shared_ptr<const GenomeReference>& genome_db_ptr) {

  std::stringstream ss;

  std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db_ptr->getContigSequence(contig_id);
  if (not contig_opt) {

    ExecEnv::log().error("outputGenomeGene(), Unexpected could not find Contig: {}", contig_id);
    return "<error>";

  }
  if (offset + region_size >= contig_opt.value()->sequence().length()) {

    ExecEnv::log().error("outputGenomeGene(), Offset: {} + Region Size: {} exceed the Contig: {} size: {}",
                         offset, region_size, contig_id, contig_opt.value()->sequence().length());
    return "<error>";

  }

  DNA5SequenceLinear mutant_sequence;
  DNA5SequenceLinear reference_sequence;
  OffsetVariantMap variant_map;

  if (not genome_variant_ptr->getSortedVariants( contig_opt.value()->contigId(),
                                                 VariantSequence::HAPLOID_PHASED,
                                                 offset,
                                                 offset+region_size,
                                                 variant_map)) {

    ExecEnv::log().warn("GenomicSequence::outputGenomeRegion, Problem retrieving variants, genome: {}, contig: {}",
                        genome_variant_ptr->genomeId(), contig_opt.value()->contigId());
    return "<error>";

  }

  if (GenomeMutation::mutantRegion( contig_id,
                                    offset,
                                    region_size,
                                    genome_db_ptr,
                                    variant_map,
                                    reference_sequence,
                                    mutant_sequence)) {

    double distance = static_cast<double>(dna_distance_metric->linear_distance(reference_sequence, mutant_sequence));

    ss << genome_variant_ptr->genomeId() << delimiter;
    ss << contig_id << delimiter;
    ss << contig_opt.value()->sequence().length() << delimiter;
    ss << offset << delimiter;
    ss << region_size << delimiter;
    ss << distance << delimiter;
    ss << SequenceComplexity::relativeCpGIslands(reference_sequence) << delimiter;
    std::vector<std::pair<DNA5::Alphabet, size_t>> symbol_count = reference_sequence.countSymbols();
    for (auto const& count : symbol_count) {

      ss << count.second << delimiter;

    }

  } else {

    ExecEnv::log().error("outputGenomeGene(), Unexpected error mutating Genome: {}, Contig: {}, Offset: {}, Size: {}",
                         genome_variant_ptr->genomeId(), contig_id, offset, region_size);
    return "<error>";

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

    ExecEnv::log().error("mutateGenomeRegion(), Genome: {} not found", genome);
    return false;

  }

}

bool kgl::GenomicSequence::mutateGenomeRegion(const ContigId_t& contig,
                                              const ContigOffset_t offset,
                                              const ContigSize_t region_size,
                                              const std::shared_ptr<const GenomeDB>& genome_variant_ptr,
                                              const std::shared_ptr<const GenomeReference>& genome_db_ptr,
                                              const std::string& fasta_file) {

  std::shared_ptr<const LinearDNASequenceDistance> dna_distance_metric(std::make_shared<const LevenshteinGlobal>());
  std::vector<WriteFastaSequence> dna_seq_vector;
  DNA5SequenceLinear mutant_sequence;
  DNA5SequenceLinear reference_sequence;
  OffsetVariantMap variant_map;

  if (not genome_variant_ptr->getSortedVariants( contig,
                                                 VariantSequence::HAPLOID_PHASED,
                                                 offset,
                                                 offset+region_size,
                                                 variant_map)) {

    ExecEnv::log().warn("GenomicSequence::mutateGenomeRegion Problem retrieving variants, genome: {}, contig: {}",
                        genome_variant_ptr->genomeId(), contig);
    return false;

  }


  if (GenomeMutation::mutantRegion( contig,
                                    offset,
                                    region_size,
                                    genome_db_ptr,
                                    variant_map,
                                    reference_sequence,
                                    mutant_sequence)) {


    DNA5SequenceCoding ref_reverse_complement = SequenceOffset::codingSequence(reference_sequence, StrandSense::REVERSE);

    CompareDistance_t score = dna_distance_metric->linear_distance(reference_sequence, mutant_sequence);

    ExecEnv::log().info("Genome: {}, Contig: {}, Offset: {} Size: {} score: {}",
                        genome_variant_ptr->genomeId(), contig, offset, region_size, score);

    auto reference_sequence_ptr = std::make_shared<DNA5SequenceLinear>(std::move(reference_sequence));
    auto mutant_sequence_ptr = std::make_shared<DNA5SequenceLinear>(std::move(mutant_sequence));
    auto ref_reverse_complement_ptr = std::make_shared<DNA5SequenceCoding>(std::move(ref_reverse_complement));

    std::stringstream ref_fasta_ss;
    ref_fasta_ss << "reference_" << contig << "_+_" << offset << "_" << region_size;
    WriteFastaSequence fasta_reference(ref_fasta_ss.str(), "", reference_sequence_ptr);
    dna_seq_vector.push_back(fasta_reference);

    ref_fasta_ss.str("");
    ref_fasta_ss << "reference_" << contig << "_-_" << offset << "_" << region_size;
    WriteFastaSequence rev_fasta_reference(ref_fasta_ss.str(), "", ref_reverse_complement_ptr);
    dna_seq_vector.push_back(rev_fasta_reference);

    std::stringstream mutant_fasta_ss;
    mutant_fasta_ss << "mutant_" << contig << "_" << offset << "_" << region_size;

    WriteFastaSequence fasta_mutant(mutant_fasta_ss.str(), "", mutant_sequence_ptr);
    dna_seq_vector.push_back(fasta_mutant);

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

