//
// Created by kellerberrin on 28/12/20.
//

#ifndef KGL_GENOME_GENOME_H
#define KGL_GENOME_GENOME_H


#include "kgl_genome_contig.h"

#include <memory>
#include <string>
#include <vector>
#include <map>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase - A map of contigs defining the genome of an organism.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeContigMap = std::map<ContigId_t, std::shared_ptr<ContigReference>>;

class GenomeReference {

public:

  explicit GenomeReference(const GenomeId_t& genome_id) : _genome_id(genome_id) {}
  GenomeReference(const GenomeReference&) = default;
  ~GenomeReference() = default;

  GenomeReference& operator=(const GenomeReference&) = default;

  // High level function creates a genome database.
  [[nodiscard]] static std::shared_ptr<GenomeReference> createGenomeDatabase(const RuntimeProperties& runtime_options,
                                                                             const GenomeId_t& organism);
  // Organism identifier
  [[nodiscard]] const GenomeId_t& genomeId() const { return _genome_id; }

  // Return false if contig already exists.
  [[nodiscard]] bool addContigSequence(const ContigId_t& contig, const std::string& description, std::shared_ptr<DNA5SequenceContig> sequence_ptr);
  // Returns false if key not found.
  [[nodiscard]] std::optional<std::shared_ptr<const ContigReference>> getContigSequence(const ContigId_t& contig) const;

  void setTranslationTable(const std::string& table);

  [[nodiscard]] const GenomeContigMap& getMap() const { return genome_sequence_map_; }

  [[nodiscard]] size_t contigCount() const { return getMap().size(); }

  [[nodiscard]] const GeneOntology& geneOntology() const { return gene_ontology_; }

  // Given a gene sequence offset with 5' start = 0 (strand adjusted), returns a strand adjusted offset within the contig.
  [[nodiscard]] bool contigOffset( const ContigId_t& contig_id,
                                   const FeatureIdent_t& gene_id,
                                   const FeatureIdent_t& sequence_id,
                                   ContigOffset_t sequence_offset,
                                   ContigOffset_t& contig_offset) const;

  // Creates a genome database object.
  // The fasta and gff files must be specified and present.
  // The gaf file is optional (empty string if omitted)
  // The translation Amino Acid table is optional (empty string if omitted).
  // Note that different translation tables can be specified for individual contigs if required.
  [[nodiscard]] static std::shared_ptr<GenomeReference> createGenomeDatabase(const GenomeId_t& organism,
                                                                             const std::string& fasta_file,
                                                                             const std::string& gff_file,
                                                                             const std::string& gaf_file,
                                                                             const std::string& translation_table);

private:

  const GenomeId_t _genome_id;
  GenomeContigMap genome_sequence_map_;
  GeneOntology gene_ontology_;

  void createVerifyGenomeDatabase();
  void createVerifyAuxillary();

  // Reads auxiliary genome information about the database. Promoter sites, motifs, tss etc.
  [[nodiscard]] bool readGenomeAuxiliary(const RuntimeProperties& runtime_options);
  // Read the auxillary genome database features.
  void readAuxillary(const std::string& tss_gff_file);

};



}   // end namespace




#endif //KGL_GENOME_GENOME_H
