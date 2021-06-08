//
// Created by kellerberrin on 28/12/20.
//

#ifndef KGL_GENOME_GENOME_H
#define KGL_GENOME_GENOME_H


#include "kgl_genome_contig.h"
#include "kgl_resource_db.h"

#include <memory>
#include <string>
#include <vector>
#include <map>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase - A map of contigs defining the genome of an organism.
// This object is passed as a requested resource into an analysis package.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeContigMap = std::map<ContigId_t, std::shared_ptr<ContigReference>>;

class GenomeReference : public ResourceBase {

public:

  explicit GenomeReference(const GenomeId_t& genome_id) : _genome_id(genome_id) {}
  GenomeReference(const GenomeReference&) = default;
  ~GenomeReference() override = default;

  // Resource type identifier.
  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::GENOME_DATABASE; }
  // Organism identifier
  [[nodiscard]] const GenomeId_t& genomeId() const { return _genome_id; }

  // Return false if contig already exists.
  [[nodiscard]] bool addContigSequence(const ContigId_t& contig, const std::string& description, std::shared_ptr<DNA5SequenceContig> sequence_ptr);
  // Returns false if key not found.
  [[nodiscard]] std::optional<std::shared_ptr<const ContigReference>> getContigSequence(const ContigId_t& contig) const;

  void setTranslationTable(const std::string& table);

  [[nodiscard]] const GenomeContigMap& getMap() const { return genome_sequence_map_; }

  [[nodiscard]] const GeneOntology& geneOntology() const { return gene_ontology_; }

  // Given a gene sequence offset with 5' start = 0 (strand adjusted), returns a strand adjusted offset within the contig.
  [[nodiscard]] bool contigOffset( const ContigId_t& contig_id,
                                   const FeatureIdent_t& gene_id,
                                   const FeatureIdent_t& sequence_id,
                                   ContigOffset_t sequence_offset,
                                   ContigOffset_t& contig_offset) const;

  // Creates a genome database object.
  // The fasta and gff files must be specified and present.
  // The gaf file and id files are optional (empty string if omitted)
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

};



}   // end namespace




#endif //KGL_GENOME_GENOME_H
