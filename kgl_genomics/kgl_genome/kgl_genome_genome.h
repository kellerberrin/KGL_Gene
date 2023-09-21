//
// Created by kellerberrin on 28/12/20.
//

#ifndef KGL_GENOME_GENOME_H
#define KGL_GENOME_GENOME_H


#include "kgl_genome_contig.h"
#include "kgl_properties_resource.h"

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

  explicit GenomeReference(const GenomeId_t& genome_id) : ResourceBase(ResourceProperties::GENOME_RESOURCE_ID_, genome_id) {}
  GenomeReference(const GenomeReference&) = default;
  ~GenomeReference() override = default;

  // Organism identifier and resource identifier
  [[nodiscard]] const GenomeId_t& genomeId() const { return resourceIdent(); }

  // ReturnType false if contig already exists.
  [[nodiscard]] bool addContigSequence(const ContigId_t& contig, const std::string& description, std::shared_ptr<DNA5SequenceLinear> sequence_ptr);
  // Returns false if key not found.
  [[nodiscard]] std::optional<std::shared_ptr<const ContigReference>> getContigSequence(const ContigId_t& contig) const;

  void setTranslationTable(const std::string& table);

  [[nodiscard]] const GenomeContigMap& getMap() const { return genome_sequence_map_; }

  [[nodiscard]] const GeneOntology& geneOntology() const { return gene_ontology_; }

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

  // Compares two genome references for equality (used for testing).
  bool equivalent(const GenomeReference& lhs) const;

private:

  GenomeContigMap genome_sequence_map_;
  GeneOntology gene_ontology_;

  void createVerifyGenomeDatabase();

};



}   // end namespace




#endif //KGL_GENOME_GENOME_H
