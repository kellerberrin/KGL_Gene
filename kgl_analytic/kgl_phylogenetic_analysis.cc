//
// Created by kellerberrin on 17/11/17.
//

#include "kgl_sequence_distance.h"
#include "kgl_sequence_compare.h"
#include "kgl_genome_aux_csv.h"
#include "kgl_upgma_node.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_sequence_offset.h"
#include "kgl_phylogenetic_gene.h"
#include "kgl_sequence_complexity.h"

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
                                             std::shared_ptr<DNA5SequenceCoding>& mutant_sequence) {

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

  std::shared_ptr<DNA5SequenceLinear> linear_mutant_sequence;
  std::shared_ptr<DNA5SequenceLinear> linear_reference_sequence;
  OffsetVariantMap variant_map;

  if (genome_variant->mutantRegion(contig_id,
                                   ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                   offset_5_prime,
                                   size_5_prime,
                                   genome_db,
                                   variant_map,
                                   linear_reference_sequence,
                                   linear_mutant_sequence)) {

    reference_sequence = SequenceOffset::codingSequence(linear_reference_sequence, coding_sequence_ptr->strand());
    mutant_sequence = SequenceOffset::codingSequence(linear_mutant_sequence, coding_sequence_ptr->strand());

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
                                             std::shared_ptr<DNA5SequenceCoding>& mutant_sequence) {

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


  std::shared_ptr<DNA5SequenceLinear> linear_mutant_sequence;
  std::shared_ptr<DNA5SequenceLinear> linear_reference_sequence;
  OffsetVariantMap variant_map;

  if (genome_variant->mutantRegion(contig_id,
                                   ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                   offset_3_prime,
                                   size_3_prime,
                                   genome_db,
                                   variant_map,
                                   linear_reference_sequence,
                                   linear_mutant_sequence)) {

    reference_sequence = SequenceOffset::codingSequence(linear_reference_sequence, coding_sequence_ptr->strand());
    mutant_sequence = SequenceOffset::codingSequence(linear_mutant_sequence, coding_sequence_ptr->strand());

  } else {

    ExecEnv::log().warn("No valid 3 prime sequence for contig: {}, offset: {}, size: {}", contig_id, offset_3_prime, size_3_prime);
    return false;

  }

  return true;

}


bool kgl::ApplicationAnalysis::outputRegionCSV(const std::string &file_name,
                                                 std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric,
                                                 std::shared_ptr<const GlobalAminoSequenceDistance> amino_distance_metric,
                                                 std::shared_ptr<const GenomeDatabase> genome_db,
                                                 std::shared_ptr<const PhasedPopulation> pop_variant_ptr) {

  const char CSV_delimiter = ',';
  // open the file.
  std::fstream out_file(file_name, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  out_file << GeneAnalysis::outputRegionHeader(CSV_delimiter) << CSV_delimiter;
  out_file << outputRegionHeader(CSV_delimiter) << '\n';

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

              ExecEnv::log().vinfo("Sequence: {} end offset: {} overlaps sequence: {} begin offset: {}",
                                  previous_seq_ptr->getCDSParent()->id(),
                                  previous_seq_ptr->end(),
                                  sequence.second->getCDSParent()->id(),
                                  sequence.second->start());

              front_porch_size = 0;

            }

          }

          out_file << GeneAnalysis::outputGenomeRegion(CSV_delimiter,
                                                       dna_distance_metric,
                                                       contig.first,
                                                       front_porch_offset,
                                                       front_porch_size,
                                                       genome_variant.second,
                                                       genome_db);
          out_file << CSV_delimiter;
          out_file << outputSequence(CSV_delimiter,
                                     dna_distance_metric,
                                     amino_distance_metric,
                                     sequence.second,
                                     genome_db,
                                     genome_variant.second);
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


bool kgl::ApplicationAnalysis::outputDNASequenceCSV(const std::string &file_name,
                                                    SequenceAnalysisType analysis_type,
                                                    std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric,
                                                    std::shared_ptr<const GenomeDatabase> genome_db,
                                                    std::shared_ptr<const PhasedPopulation> pop_variant_ptr) {

  const char CSV_delimiter = ',';
  // open the file.
  std::ofstream out_file(file_name);
  if (not out_file.good()) {

    ExecEnv::log().error("outputDNASequenceCSV, Cannot open output CSV file: {}", file_name);
    return false;

  }

  out_file << outputSequenceHeader(CSV_delimiter, pop_variant_ptr) << '\n';

  for (auto contig : genome_db->getMap()) {

    ExecEnv::log().info("outputDNASequenceCSV(), Processing contig: {}", contig.first);

    size_t sequence_count = 0;

    for (auto gene : contig.second->getGeneMap()) {

      const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene.second);
      for (auto sequence : coding_seq_ptr->getMap()) {


        out_file << contig.second->contigSize() << CSV_delimiter;
        out_file << gene.second->id() << CSV_delimiter;
        out_file << sequence.first << CSV_delimiter;
        out_file << sequence.second->start() << CSV_delimiter;
        out_file << sequence.second->codingNucleotides() << CSV_delimiter;

        // Valid ORF
        std::shared_ptr<DNA5SequenceCoding> coding_dna_sequence;
        if (contig.second->getDNA5SequenceCoding(sequence.second, coding_dna_sequence)) {

          out_file << (contig.second->verifyDNACodingSequence(coding_dna_sequence) ? "1" : "0") << CSV_delimiter;

        } else {

          out_file << "0" << CSV_delimiter;

        }

        out_file << SequenceComplexity::alphabetEntropy<CodingDNA5>(coding_dna_sequence) << CSV_delimiter;
        out_file << SequenceComplexity::complexityLempelZiv(DNA5SequenceLinear::downConvertToLinear(coding_dna_sequence)) << CSV_delimiter;
        out_file << SequenceComplexity::relativeCpGIslands(DNA5SequenceLinear::downConvertToLinear(coding_dna_sequence)) << CSV_delimiter;  // GC count.
        double A_prop;
        double C_prop;
        double G_prop;
        double T_prop;
        SequenceComplexity::proportionNucleotides(DNA5SequenceLinear::downConvertToLinear(coding_dna_sequence), A_prop, C_prop, G_prop, T_prop);
        out_file << A_prop << CSV_delimiter;
        out_file << C_prop << CSV_delimiter;
        out_file << G_prop << CSV_delimiter;
        out_file << T_prop << CSV_delimiter;

        std::shared_ptr<const OntologyRecord> gene_ontology_ptr;
        if (genome_db->geneOntology().getGafFeatureVector(gene.second->id(), gene_ontology_ptr)) {

          out_file << gene_ontology_ptr->symbolicReference() << CSV_delimiter;
          out_file << gene_ontology_ptr->altSymbolicReference() << CSV_delimiter;
          out_file << "\"" << gene_ontology_ptr->description() << "\"" << CSV_delimiter;

          std::set<std::string> uniqueGOrecords;
          for (auto gorec : gene_ontology_ptr->goRecords()) {

            uniqueGOrecords.insert(gorec.first);

          }
          std::string GOrecords;
          for (auto gorec : uniqueGOrecords) {

            GOrecords += gorec;
            GOrecords += ";";

          }

          out_file << "\"" << GOrecords << "\"" << CSV_delimiter;

        } else {

          out_file << CSV_delimiter;
          out_file << CSV_delimiter;
          std::vector<std::string> description_vec;
          gene.second->getAttributes().getDescription(description_vec);
          std::string description;
          for (const auto &desc : description_vec) {

            description += desc;
            description += ";";

          }

          out_file << "\"" << description << "\"" << CSV_delimiter;
          out_file << CSV_delimiter;

        }


        for( auto genome_variant : pop_variant_ptr->getMap()) {

          std::shared_ptr<DNA5SequenceCoding> reference_sequence;
          std::shared_ptr<DNA5SequenceCoding> mutant_sequence;
          OffsetVariantMap variant_map;

          if (genome_variant.second->mutantCodingDNA( contig.first,
                                                      ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                                      gene.second->id(),
                                                      sequence.first,
                                                      genome_db,
                                                      variant_map,
                                                      reference_sequence,
                                                      mutant_sequence)) {


            switch(analysis_type) {

              case SequenceAnalysisType::DNA: {

                CompareDistance_t DNA_distance;
                DNA_distance = dna_distance_metric->distance(reference_sequence, mutant_sequence);
                out_file << DNA_distance << CSV_delimiter;

              }
                break;

              case SequenceAnalysisType::SIZE:
                out_file << mutant_sequence->length() << CSV_delimiter;
                break;

              case SequenceAnalysisType::VARIANT:
                out_file << variant_map.size() << CSV_delimiter;
                break;

              case SequenceAnalysisType::SNP: {

                size_t snp_count = 0;
                for (auto variant : variant_map) {

                  if (variant.second->isSNP()) ++snp_count;

                }

                out_file << variant_map.size() << CSV_delimiter;

              }
                break;

              case SequenceAnalysisType::ENTROPY:
                out_file << SequenceComplexity::alphabetEntropy<CodingDNA5>(mutant_sequence) << CSV_delimiter;
                break;

              case SequenceAnalysisType::LEMPEL_ZIV:
                out_file << SequenceComplexity::complexityLempelZiv(
                DNA5SequenceLinear::downConvertToLinear(mutant_sequence)) << CSV_delimiter;
                break;

            }


          } else {

            ExecEnv::log().error("outputDNASequenceCSV(), Error Processing sequence: {}", sequence.first);
            return false;

          }

          ++sequence_count;

        }

        out_file << '\n';

      }

    }

    ExecEnv::log().info("outputDNASequenceCSV(), mutated: {} sequences.", sequence_count);

  }

  return out_file.good();

}


bool kgl::ApplicationAnalysis::outputAminoSequenceCSV(const std::string &file_name,
                                                      std::shared_ptr<const GlobalAminoSequenceDistance> amino_distance_metric,
                                                      std::shared_ptr<const GenomeDatabase> genome_db,
                                                      std::shared_ptr<const PhasedPopulation> pop_variant_ptr) {

  const char CSV_delimiter = ',';
  // open the file.
  std::ofstream out_file(file_name);
  if (not out_file.good()) {

    ExecEnv::log().error("outputAminoSequenceCSV, Cannot open output CSV file: {}", file_name);
    return false;

  }

  out_file << outputSequenceHeader(CSV_delimiter, pop_variant_ptr) << '\n';

  for (auto contig : genome_db->getMap()) {

    ExecEnv::log().info("outputAminoSequenceCSV(), Processing contig: {}", contig.first);

    size_t sequence_count = 0;

    for (auto gene : contig.second->getGeneMap()) {

      const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene.second);
      for (auto sequence : coding_seq_ptr->getMap()) {


        out_file << contig.second->contigSize() << CSV_delimiter;
        out_file << gene.second->id() << CSV_delimiter;
        out_file << sequence.first << CSV_delimiter;
        out_file << sequence.second->start() << CSV_delimiter;
        out_file << sequence.second->codingNucleotides() << CSV_delimiter;

        std::shared_ptr<DNA5SequenceCoding> coding_dna_sequence;
        if (contig.second->getDNA5SequenceCoding(sequence.second, coding_dna_sequence)) {

          out_file << (contig.second->verifyDNACodingSequence(coding_dna_sequence) ? "1" : "0") << CSV_delimiter;

        } else {

          out_file << "0" << CSV_delimiter;

        }

        std::shared_ptr<AminoSequence> amino_reference  = contig.second->getAminoSequence(coding_dna_sequence);
        out_file << SequenceComplexity::alphabetEntropy<AminoAcid>(amino_reference) << CSV_delimiter;
        out_file << SequenceComplexity::complexityLempelZiv(DNA5SequenceLinear::downConvertToLinear(coding_dna_sequence)) << CSV_delimiter;
        out_file << SequenceComplexity::relativeCpGIslands(DNA5SequenceLinear::downConvertToLinear(coding_dna_sequence)) << CSV_delimiter;  // GC count.
        double A_prop;
        double C_prop;
        double G_prop;
        double T_prop;
        SequenceComplexity::proportionNucleotides(DNA5SequenceLinear::downConvertToLinear(coding_dna_sequence), A_prop, C_prop, G_prop, T_prop);
        out_file << A_prop << CSV_delimiter;
        out_file << C_prop << CSV_delimiter;
        out_file << G_prop << CSV_delimiter;
        out_file << T_prop << CSV_delimiter;

        std::shared_ptr<const OntologyRecord> gene_ontology_ptr;
        if (genome_db->geneOntology().getGafFeatureVector(gene.second->id(), gene_ontology_ptr)) {

          out_file << gene_ontology_ptr->symbolicReference() << CSV_delimiter;
          out_file << gene_ontology_ptr->altSymbolicReference() << CSV_delimiter;
          out_file << "\"" << gene_ontology_ptr->description() << "\"" << CSV_delimiter;

          std::set<std::string> uniqueGOrecords;
          for (auto gorec : gene_ontology_ptr->goRecords()) {

            uniqueGOrecords.insert(gorec.first);

          }
          std::string GOrecords;
          for (auto gorec : uniqueGOrecords) {

            GOrecords += gorec;
            GOrecords += ";";

          }

          out_file << "\"" << GOrecords << "\"" << CSV_delimiter;

        } else {

          out_file << CSV_delimiter;
          out_file << CSV_delimiter;
          std::vector<std::string> description_vec;
          gene.second->getAttributes().getDescription(description_vec);
          std::string description;
          for (const auto &desc : description_vec) {

            description += desc;
            description += ";";

          }

          out_file << "\"" << description << "\"" << CSV_delimiter;
          out_file << CSV_delimiter;

        }


        for( auto genome_variant : pop_variant_ptr->getMap()) {

          OffsetVariantMap variant_map;
          std::shared_ptr<AminoSequence> amino_reference_seq;
          std::shared_ptr<AminoSequence> amino_mutant;
          if (genome_variant.second->mutantProteins(contig.first,
                                                    ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                                    gene.second->id(),
                                                    sequence.first,
                                                    genome_db,
                                                    variant_map,
                                                    amino_reference_seq,
                                                    amino_mutant)) {

            switch (contig.second->proteinSequenceAnalysis(amino_mutant)) {

              case ProteinSequenceAnalysis::VALID_SEQUENCE: {

                CompareDistance_t amino_distance = amino_distance_metric->distance(amino_reference_seq, amino_mutant);
                out_file << amino_distance << CSV_delimiter;

              }
                break;

              case ProteinSequenceAnalysis::NO_START_CODON:
                out_file << -1 << CSV_delimiter;
                break;

              case ProteinSequenceAnalysis::NONSENSE_MUTATION: {

                size_t valid_seq_size = contig.second->proteinSequenceSize(amino_mutant) + 1;  // +1 to include the stop codon.
                double proportion = (static_cast<double>(valid_seq_size) * -100.0) / static_cast<double>(amino_reference_seq->length());
                out_file << proportion << CSV_delimiter;

              }
                break;

              case ProteinSequenceAnalysis::NO_STOP_CODON:
                out_file << -3 << CSV_delimiter;
                break;

            }

          } else {

            ExecEnv::log().error("outputSequenceCSV(), Error Processing sequence: {}", sequence.first);
            return false;

          }

          ++sequence_count;

        }

        out_file << '\n';

      }

    }

    ExecEnv::log().info("outputAminoSequenceCSV(), mutated: {} sequences.", sequence_count);

  }

  return out_file.good();

}



std::string kgl::ApplicationAnalysis::outputSequenceHeader(char delimiter,
                                                           std::shared_ptr<const PhasedPopulation> pop_variant_ptr) {

  std::stringstream ss;

  ss << "Contig" << delimiter;
  ss << "Gene" << delimiter;
  ss << "Sequence" << delimiter;
  ss << "Offset" << delimiter;
  ss << "DNALength" << delimiter;
  ss << "ValidORF" << delimiter;
  ss << "Entropy" << delimiter;
  ss << "LempelZiv" << delimiter;
  ss << "CpG" << delimiter;
  ss << "A_prop" << delimiter;
  ss << "C_prop" << delimiter;
  ss << "G_prop" << delimiter;
  ss << "T_prop" << delimiter;
  ss << "Symbolic" << delimiter;
  ss << "AltSymbolic" << delimiter;
  ss << "Description" << delimiter;
  ss << "GO_Records" << delimiter;


  for (auto genome : pop_variant_ptr->getMap()) {

    ss << genome.first << delimiter;

  }

  return ss.str();

}


bool kgl::ApplicationAnalysis::outputAminoMutationCSV(const std::string &file_name,
                                                 const ContigId_t& contig_id,
                                                 const FeatureIdent_t& gene_id,
                                                 const FeatureIdent_t& sequence_id,
                                                 std::shared_ptr<const GenomeDatabase> genome_db,
                                                 std::shared_ptr<const PhasedPopulation> pop_variant_ptr) {

  const char CSV_delimiter = ',';
  // open the file.
  std::fstream out_file(file_name, std::fstream::out);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  for( auto genome_variant : pop_variant_ptr->getMap()) {

    ExecEnv::log().info("outputMutationCSV(), Processing genome: {}", genome_variant.first);
    size_t sequence_count = 0;

    sequence_count++;
    OffsetVariantMap variant_map;
    std::shared_ptr<AminoSequence> amino_reference_seq;
    std::shared_ptr<AminoSequence> amino_mutant;
    if (genome_variant.second->mutantProteins(contig_id,
                                              ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                              gene_id,
                                              sequence_id,
                                              genome_db,
                                              variant_map,
                                              amino_reference_seq,
                                              amino_mutant)) {

        EditVector edit_vector;
        SequenceComparison().editDNAItems(amino_reference_seq->getSequenceAsString(),
                                          amino_mutant->getSequenceAsString(),
                                          edit_vector);

      for (auto edit_item : edit_vector) {

        std::stringstream ss;
        ss << genome_variant.first << CSV_delimiter;
        ss << gene_id << CSV_delimiter;
        ss << sequence_id << CSV_delimiter;
        ss << edit_item.reference_char << CSV_delimiter;
        ss << offsetOutput(edit_item.reference_offset, VariantOutputIndex::START_1_BASED) << CSV_delimiter;
        ss << edit_item.mutant_char << '\n';
        out_file << ss.str();

      } // for mutant

    } // if mutation

    ExecEnv::log().info("outputMutantCSV(), Genome: {} mutated: {} sequences.", genome_variant.first, sequence_count);

  } // for genome

  return out_file.good();

}



bool kgl::ApplicationAnalysis::outputDNAMutationCSV(const std::string &file_name,
                                                    const ContigId_t& contig_id,
                                                    const FeatureIdent_t& gene_id,
                                                    const FeatureIdent_t& sequence_id,
                                                    std::shared_ptr<const GenomeDatabase> genome_db,
                                                    std::shared_ptr<const PhasedPopulation> pop_variant_ptr,
                                                    const GenomeAuxData& aux_Pf3k_data) {

  const char CSV_delimiter = ',';
  // open the file.
  std::fstream out_file(file_name, std::fstream::out);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file: {}", file_name);
    return false;

  }

  std::fstream variant_file(file_name + ".variant", std::fstream::out);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output Variant file: {}", file_name + ".variant");
    return false;

  }

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().error("outputMutationCSV(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  using SNPMap = std::map<std::string, MutationItem>;
  SNPMap master_SNP_List;
  struct GenomeMap {

    std::string genome;
    std::string location_date;
    SNPMap snp_map;

  };

  std::vector<GenomeMap> genome_vector;
  for( auto genome_variant : pop_variant_ptr->getMap()) {

    ExecEnv::log().info("outputMutationCSV(), Processing genome: {}", genome_variant.first);
    size_t sequence_count = 0;

    sequence_count++;
    OffsetVariantMap variant_map;
    std::shared_ptr<DNA5SequenceCoding> reference_sequence;
    std::shared_ptr<DNA5SequenceCoding> mutant_sequence;
    if (genome_variant.second->mutantCodingDNA( contig_id,
                                                ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                                gene_id,
                                                sequence_id,
                                                genome_db,
                                                variant_map,
                                                reference_sequence,
                                                mutant_sequence)) {

      for (auto variant : variant_map) {

        variant_file << variant.second->output(CSV_delimiter, VariantOutputIndex::START_1_BASED, false);

      }

      EditVector edit_vector;
      SequenceComparison().editDNAItems(reference_sequence->getSequenceAsString(),
                                        mutant_sequence->getSequenceAsString(),
                                        edit_vector);

      MutationEditVector mutation_edit_vector;
      MutationItem mutation_item;
      for (auto edit_item : edit_vector) {

        ContigOffset_t codon_index = static_cast<size_t>(edit_item.reference_offset / 3);

        std::shared_ptr<const Codon> mutant_codon(std::make_shared<Codon>(mutant_sequence, codon_index));
        std::shared_ptr<const Codon> ref_codon(std::make_shared<Codon>(reference_sequence, codon_index));

        mutation_item.DNA_mutation = edit_item;
        mutation_item.reference_codon = ref_codon->getSequenceAsString();
        mutation_item.mutation_codon = mutant_codon->getSequenceAsString();

        mutation_item.amino_mutation.reference_char = AminoAcid::convertToChar(contig_ptr->getAminoAcid(*ref_codon));
        mutation_item.amino_mutation.reference_offset = codon_index;
        mutation_item.amino_mutation.mutant_char = AminoAcid::convertToChar(contig_ptr->getAminoAcid(*mutant_codon));
        mutation_item.contig_id = contig_ptr->contigId();
        ContigOffset_t contig_offset;
        if (not genome_db->contigOffset(contig_id, gene_id, sequence_id, mutation_item.DNA_mutation.reference_offset, contig_offset)) {

          ExecEnv::log().error("Edit Item {}{}{} sequence offset out of range",
                               mutation_item.DNA_mutation.reference_char, mutation_item.DNA_mutation.reference_offset, mutation_item.DNA_mutation.mutant_char);

        }
        mutation_item.contig_offset = contig_offset;
        mutation_edit_vector.mutation_vector.push_back(mutation_item);

      }

      if (not mutation_edit_vector.hasIndel()) { // Only SNPs.

        GenomeMap genome_map;

        genome_map.genome = genome_variant.first;
        if (aux_Pf3k_data.isFieldSample(genome_variant.first)) {

          genome_map.location_date = aux_Pf3k_data.locationDate(genome_variant.first);

        } else {

          genome_map.location_date = genome_variant.first;

        }

        for (auto edit_item : mutation_edit_vector.mutation_vector) {

          master_SNP_List[edit_item.mapKey()] = edit_item;
          genome_map.snp_map[edit_item.mapKey()] = edit_item;
          ExecEnv::log().info("Edit Item Contig: {}, Contig Offset:{}, DNA Offset:{}",
                              edit_item.contig_id, edit_item.contig_offset, master_SNP_List[edit_item.mapKey()].DNA_mutation.reference_offset);

        }

        genome_vector.push_back(genome_map);

      }

    } // if mutation

    ExecEnv::log().info("outputMutantCSV(), Genome: {} mutated: {} sequences.", genome_variant.first, sequence_count);

  } // for genome

// Write file header.

  out_file << "Genome" << CSV_delimiter << "LocationDate" << CSV_delimiter;

  for (auto DNA_Item : master_SNP_List) {

    out_file << DNA_Item.second.contig_id << " "
             << offsetOutput(DNA_Item.second.contig_offset, VariantOutputIndex::START_1_BASED) << " "
             << DNA_Item.second.DNA_mutation.reference_char
             << offsetOutput(DNA_Item.second.DNA_mutation.reference_offset, VariantOutputIndex::START_1_BASED)
             << DNA_Item.second.DNA_mutation.mutant_char
             << " " << DNA_Item.second.reference_codon << "-" << DNA_Item.second.mutation_codon
             << " " << DNA_Item.second.amino_mutation.reference_char
             << offsetOutput(DNA_Item.second.amino_mutation.reference_offset, VariantOutputIndex::START_1_BASED)
             << DNA_Item.second.amino_mutation.mutant_char << CSV_delimiter;

  }
  out_file << std::endl;

// Write data.

  for (auto genome_item : genome_vector) {

    out_file << genome_item.genome << CSV_delimiter << genome_item.location_date << CSV_delimiter;

    for (auto DNA_Item : master_SNP_List) {

      auto result = genome_item.snp_map.find(DNA_Item.first);
      if (result == genome_item.snp_map.end()) {

        out_file << 0 << CSV_delimiter;

      } else {

        out_file << 1 << CSV_delimiter;

      }

    }

    out_file << std::endl;

  }

  return out_file.good();

}



std::string kgl::ApplicationAnalysis::outputRegionHeader(char delimiter) {

  std::stringstream ss;

  ss << "Genome" << delimiter;
  ss << "Contig" << delimiter;
  ss << "ContigLength" << delimiter;
  ss << "Sequence" << delimiter;
  ss << "ContigOffset" << delimiter;
  ss << "DNASize" << delimiter;
  ss << "CpG" << delimiter;
  ss << "A_prop" << delimiter;
  ss << "C_prop" << delimiter;
  ss << "G_prop" << delimiter;
  ss << "T_prop" << delimiter;
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
                                                     std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric,
                                                     std::shared_ptr<const GlobalAminoSequenceDistance> amino_distance_metric,
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
  double average_DNA_score = 0;
  bool valid_reference = false;
  OffsetVariantMap variant_map;

  std::shared_ptr<DNA5SequenceCoding> reference_sequence;
  std::shared_ptr<DNA5SequenceCoding> mutant_sequence;
  if (genome_variant->mutantCodingDNA( contig_ptr->contigId(),
                                       ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                       gene_id,
                                       sequence_id,
                                       genome_db,
                                       variant_map,
                                       reference_sequence,
                                       mutant_sequence)) {

    CompareDistance_t DNA_distance;
    DNA_distance = dna_distance_metric->distance(reference_sequence, mutant_sequence);
    average_DNA_score += static_cast<double>(DNA_distance);

  }

  std::shared_ptr<AminoSequence> amino_reference_seq;
  std::shared_ptr<AminoSequence> amino_mutant;
  if (genome_variant->mutantProteins(contig_ptr->contigId(),
                                     ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                     gene_id,
                                     sequence_id,
                                     genome_db,
                                     variant_map,
                                     amino_reference_seq,
                                     amino_mutant)) {

    error_flag = false;
    valid_reference = contig_ptr->verifyProteinSequence(amino_reference_seq);

    CompareDistance_t amino_distance = amino_distance_metric->distance(amino_reference_seq, amino_mutant);
    if (contig_ptr->verifyProteinSequence(amino_mutant)) {

      average_score += static_cast<double>(amino_distance);
      ++valid_paths;

    }
    ++mutant_paths;

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
  ss << sequence_id << delimiter;
  ss << contig_ptr->sequence().length() << delimiter;
  ss << sequence_offset << delimiter;
  ss << coding_sequence->codingNucleotides() << delimiter;
  ss << SequenceComplexity::relativeCpGIslands(DNA5SequenceLinear::downConvertToLinear(reference_sequence)) << delimiter;  // GC count.
  double A_prop;
  double C_prop;
  double G_prop;
  double T_prop;
  SequenceComplexity::proportionNucleotides(DNA5SequenceLinear::downConvertToLinear(reference_sequence), A_prop, C_prop, G_prop, T_prop);
  ss << A_prop << delimiter;
  ss << C_prop << delimiter;
  ss << G_prop << delimiter;
  ss << T_prop << delimiter;
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

