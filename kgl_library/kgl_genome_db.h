//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_GENOME_DB_H
#define KGL_GENOME_DB_H


#include <memory>
#include <string>
#include <vector>
#include <map>
#include <kgl_properties.h>
#include "kgl_genome_types.h"
#include "kgl_sequence_amino.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_contig_feature.h"
#include "kgl_genome_contig_aux.h"
#include "kgl_gaf_parser.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigFeatures - A contiguous region, the associated sequence, and all features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class ProteinSequenceAnalysis { VALID_SEQUENCE, NO_START_CODON, NONSENSE_MUTATION, NO_STOP_CODON };


class ContigFeatures {

public:

  ContigFeatures(const ContigId_t& contig_id,
                 std::shared_ptr<DNA5SequenceContig> sequence_ptr) : contig_id_(contig_id),
                                                                     sequence_ptr_(sequence_ptr) {}
  ContigFeatures(const ContigFeatures&) = default;
  ~ContigFeatures() = default;

  ContigFeatures& operator=(const ContigFeatures&) = default;

  // Add parsed features to the different feature structures.
  bool addGeneExonFeature(std::shared_ptr<Feature>& feature_ptr);
  bool addAuxFeature(std::shared_ptr<Feature>& feature_ptr);

  // false if not found.
  bool findFeatureId(const FeatureIdent_t& feature_id,
                     std::vector<std::shared_ptr<const Feature>>& feature_ptr_vec) const {
    return gene_exon_features_.findFeatureId(feature_id, feature_ptr_vec);
  }
  // false if offset is not in a gene, else (true) returns a vector of ptrs to the genes.
  bool findGenes(ContigOffset_t offset, GeneVector &gene_ptr_vec) const;
 
  const GeneMap& getGeneMap() const { return gene_exon_features_.geneMap(); }

  // Return all Aux genome features in this contig.
  const AuxContigFeatures& getAuxContigFeatures() const { return aux_contig_features_; }

  bool setTranslationTable(const std::string& table_name) { return coding_table_.settranslationTable(table_name); }
 
  std::string translationTableName() const { return coding_table_.translationTableName(); }

  const ContigId_t& contigId() const { return contig_id_; }
  const DNA5SequenceContig& sequence() const { return *sequence_ptr_; }
  std::shared_ptr<const DNA5SequenceContig> sequence_ptr() const { return sequence_ptr_; }
  ContigSize_t contigSize() const { return sequence_ptr_->length(); }


  bool verifyDNACodingSequence(std::shared_ptr<const DNA5SequenceCoding> coding_dna_ptr) const;
  bool verifyProteinSequence(std::shared_ptr<const AminoSequence> amino_sequence_ptr) const;
  ProteinSequenceAnalysis proteinSequenceAnalysis(std::shared_ptr<const AminoSequence> amino_sequence_ptr) const;
  // Returns the protein sequence as the distance in amino acids between the start codon and stop codon.
  // No start and stop codon returns 0.
  size_t proteinSequenceSize(std::shared_ptr<const AminoSequence> amino_sequence_ptr) const;

  // Given a gene id and an mRNA (sequence id) return the CDS coding sequence.
  bool getCodingSequence(const FeatureIdent_t& gene_id,
                         const FeatureIdent_t& sequence_id,
                         std::shared_ptr<const CodingSequence>& coding_sequence_ptr) const;

  // Given a CDS coding sequence, return the corresponding DNA base sequence (strand adjusted).
  bool getDNA5SequenceCoding(const std::shared_ptr<const CodingSequence>& coding_sequence_ptr,
                             std::shared_ptr<DNA5SequenceCoding>& coding_dna_ptr) const;

  // Generate Amino acid sequences using the table specified for this contig.
  std::shared_ptr<AminoSequence> getAminoSequence(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const;
  AminoAcid::Alphabet getAminoAcid(const Codon& codon) const { return coding_table_.getAmino(codon); }

  // Wire-up the contig features
  void verifyFeatureHierarchy();
  void verifyAuxillaryHierarchy();
  void verifyCDSPhasePeptide();

private:

  ContigId_t contig_id_;
  std::shared_ptr<DNA5SequenceContig> sequence_ptr_;  // The contig unstranded DNA sequence.
  GeneExonFeatures gene_exon_features_;  // All the genes and exons defined for this contig.
  AuxContigFeatures aux_contig_features_;
  TranslateToAmino coding_table_;  // Amino Acid translation table, unique for contig (e.g. mitochondria)

  // Check all gene coding sequences for start and end codons and nonsense (intermediate stop codon) mutations.
  bool verifyCodingSequences(const std::shared_ptr<const GeneFeature> gene_ptr,
                             std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) const;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase - A map of contigs defining the genome of an organism.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeContigMap = std::map<ContigId_t, std::shared_ptr<ContigFeatures>>;

class GenomeDatabase {

public:

  explicit GenomeDatabase(const GenomeId_t& genome_id) : _genome_id(genome_id) {}
  GenomeDatabase(const GenomeDatabase&) = default;
  ~GenomeDatabase() = default;

  GenomeDatabase& operator=(const GenomeDatabase&) = default;

  // High level function creates a genome database.
  static std::shared_ptr<GenomeDatabase> createGenomeDatabase(const RuntimeProperties& runtime_options,
                                                              const GenomeId_t& organism);
  // Organism identifier
  const GenomeId_t& genomeId() const { return _genome_id; }

  // Return false if contig already exists.
  bool addContigSequence(const ContigId_t& contig, std::shared_ptr<DNA5SequenceContig> sequence_ptr);
  // Returns false if key not found.
  bool getContigSequence(const ContigId_t& contig, std::shared_ptr<const ContigFeatures>& contig_ptr) const;

  void setTranslationTable(const std::string& table);

  const GenomeContigMap& getMap() const { return genome_sequence_map_; }

  size_t contigCount() const { return getMap().size(); }

  const GeneOntology& geneOntology() const { return gene_ontology_; }

  // Given a gene sequence offset with 5' start = 0 (strand adjusted), returns a strand adjusted offset within the contig.
  bool contigOffset( const ContigId_t& contig_id,
                     const FeatureIdent_t& gene_id,
                     const FeatureIdent_t& sequence_id,
                     ContigOffset_t sequence_offset,
                     ContigOffset_t& contig_offset) const;

private:

  const GenomeId_t _genome_id;
  GenomeContigMap genome_sequence_map_;
  GeneOntology gene_ontology_;

  void createVerifyGenomeDatabase();
  void createVerifyAuxillary();
  // Creates a genome database object.
  // The fasta and gff files must be specified and present.
  // The gaf file is optional (empty string if omitted)
  // The translation Amino Acid table is optional (empty string if omitted).
  // Note that different translation tables can be specified for individual contigs.
  static std::shared_ptr<GenomeDatabase> createGenomeDatabase(const GenomeId_t& organism,
                                                              const std::string& fasta_file,
                                                              const std::string& gff_file,
                                                              const std::string& gaf_file,
                                                              const std::string& translation_table);

  // Reads auxiliary genome information about the database. Promoter sites, motifs, tss etc.
  bool readGenomeAuxiliary(const RuntimeProperties& runtime_options);
  // Read the auxillary genome database features.
  void readAuxillary(const std::string& tss_gff_file);

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeCollection - A map of different organism genomes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeMap = std::map<GenomeId_t, std::shared_ptr<const GenomeDatabase>>;

class GenomeCollection {

public:

  explicit GenomeCollection() = default;
  GenomeCollection(const GenomeCollection&) = default;
  virtual ~GenomeCollection() = default;

  GenomeCollection& operator=(const GenomeCollection&) = default;

  // High level function creates a collection of genomes.
  static std::shared_ptr<GenomeCollection> createGenomeCollection(const RuntimeProperties& runtime_options);

  // Returns false if the genome does not exist.
  bool getGenome(const GenomeId_t& genome_id, std::shared_ptr<const GenomeDatabase>& genome_variant) const;

  // Convenience Routines. Get a particular XML tagged organism (see tags below).
  std::shared_ptr<const GenomeDatabase> get3D7Genome() const { return getGenome(Pf3D7_ID_); }
  std::shared_ptr<const GenomeDatabase> getHB3Genome() const { return getGenome(PfHB3_ID_); }
  std::shared_ptr<const GenomeDatabase> getITGenome() const { return getGenome(PfIT_ID_); }
  std::shared_ptr<const GenomeDatabase> getDD2Genome() const { return getGenome(PfDD2_ID_); }
  std::shared_ptr<const GenomeDatabase> get7G8Genome() const { return getGenome(Pf7G8_ID_); }

  const GenomeMap& getMap() const { return genome_map_; }

private:

  // A map of all active genome databases.
  GenomeMap genome_map_;

  // Must match the relevant XML tag in the runtime options file.
  static constexpr const char* Pf3D7_ID_ = "Pf3D7_46";
  static constexpr const char* PfHB3_ID_ = "PfHB3_46";
  static constexpr const char* PfIT_ID_ = "PfIT_46";
  static constexpr const char* PfDD2_ID_ = "PfDd2_46";
  static constexpr const char* Pf7G8_ID_ = "Pf7G8_46";
  // Returns false if the genome already exists.
  bool addGenome(std::shared_ptr<const GenomeDatabase> genome_database);
  std::shared_ptr<const GenomeDatabase> getGenome(const std::string& GenomeID) const;

};


}   // end namespace


#endif //KGL_GENOME_DB_H
