//
// Created by kellerberrin on 3/01/18.
//

#ifndef KGL_VARIANT_DB_GENOME_H
#define KGL_VARIANT_DB_GENOME_H



#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant.h"
#include "kgl_variant_db_contig.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The variant contig map
using GenomeVariantMap = std::map<ContigId_t, std::shared_ptr<ContigVariant>>;

class GenomeVariant {

public:

  explicit GenomeVariant(const GenomeId_t& genome_id) : genome_id_(genome_id) {}
  GenomeVariant(const GenomeVariant&) = delete; // Use deep copy.
  ~GenomeVariant() = default;

  GenomeVariant& operator=(const GenomeVariant&) = delete; // Use deep copy.

  // Always use deep copy when modifying this object (filter and set operations).
  std::shared_ptr<GenomeVariant> deepCopy() const;

  const GenomeId_t& genomeId() const { return genome_id_; }
  void genomeId(const GenomeId_t& contig_id) { genome_id_ = contig_id; }

  size_t contigCount() const { return genome_variant_map_.size(); }

  bool addContigVariant(std::shared_ptr<ContigVariant>& contig_variant);
  bool getContigVariant(const ContigId_t& contig_id, std::shared_ptr<ContigVariant>& contig_variant) const;

  bool addVariant(std::shared_ptr<const Variant> variant);

  size_t size() const;

  std::shared_ptr<GenomeVariant> filterVariants(const VariantFilter& filter) const;

  const GenomeVariantMap& getMap() const { return genome_variant_map_; }

  bool isElement(const Variant& variant) const;
  std::shared_ptr<GenomeVariant> Union(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Intersection(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Difference(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;

  static std::shared_ptr<GenomeVariant> emptyGenomeVariant(const GenomeId_t& genome_id,
                                                           const std::shared_ptr<const GenomeDatabase>& genome_db);

  std::string output(char field_delimiter, VariantOutputIndex output_index, bool detail) const;
  bool outputCSV(const std::string& file_name, VariantOutputIndex output_index, bool detail) const;

  void getVariants(std::vector<std::shared_ptr<const Variant>>& variant_vector) const;

  // All contig_id variants use the zero-based half-open convention [start, end).
  // End points past the last variant; end = (last + 1).
  bool getSortedVariants(ContigId_t contig_id,
                         ContigOffset_t start,
                         ContigOffset_t end,
                         OffsetVariantMap& variant_map) const;

  // Convenience routine only returns coding variants.
  bool getCodingSortedVariants(std::shared_ptr<const CodingSequence> coding_sequence_ptr,
                               OffsetVariantMap& variant_map,
                               bool& frame_shift) const;

  const Attributes& attributes() const { return attributes_; }
  Attributes& attributes() { return attributes_; }
  void attributes(const Attributes& attributes) { attributes_ = attributes; }

  // Returns a maximum of MUTATION_HARD_LIMIT_ alternative protein mutations.
  bool mutantProteins( const ContigId_t& contig_id,
                       const FeatureIdent_t& gene_id,
                       const FeatureIdent_t& sequence_id,
                       const std::shared_ptr<const GenomeDatabase>& genome_db,
                       bool& frame_shift_mutation,
                       std::shared_ptr<AminoSequence>& reference_sequence,
                       std::vector<std::shared_ptr<AminoSequence>>& mutant_sequence_vector) const;

  // Returns a maximum of MUTATION_HARD_LIMIT_ alternative protein mutations.
  bool mutantCodingDNA( const ContigId_t& contig_id,
                        const FeatureIdent_t& gene_id,
                        const FeatureIdent_t& sequence_id,
                        const std::shared_ptr<const GenomeDatabase>& genome_db,
                        bool& frame_shift_flag,
                        std::shared_ptr<DNA5SequenceCoding>& reference_sequence,
                        std::vector<std::shared_ptr<DNA5SequenceCoding>>& mutant_sequence_vector) const;

  // Returns a maximum of MUTATION_HARD_LIMIT_ alternative mutations.
  bool mutantRegion( const ContigId_t& contig_id,
                     ContigOffset_t region_offset,
                     ContigSize_t region_size,
                     const std::shared_ptr<const GenomeDatabase>& genome_db,
                     std::shared_ptr<DNA5SequenceLinear>& reference_sequence,
                     std::vector<std::shared_ptr<DNA5SequenceLinear>>& mutant_sequence_vector) const;



private:

  GenomeId_t genome_id_;
  GenomeVariantMap genome_variant_map_;
  Attributes attributes_;

};


}   // namespace genome
}   // namespace kellerberrin


// Not in kgl:: namespace.
std::ostream & operator<<(std::ostream &os, const kellerberrin::genome::GenomeVariant& genome_variant);
std::ostream & operator<<(std::ostream &os, std::shared_ptr<const kellerberrin::genome::GenomeVariant> genome_variant_ptr);
std::ostream & operator<<(std::ostream &os, std::shared_ptr<kellerberrin::genome::GenomeVariant> genome_variant_ptr);


#endif //KGL_VARIANT_DB_GENOME_H
