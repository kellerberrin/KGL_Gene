//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_VARIANT_DB_H
#define KGL_VARIANT_DB_H


#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetVariantMap = std::multimap<ContigOffset_t, std::shared_ptr<const Variant>>;
class ContigVariant {

public:

  explicit ContigVariant(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  ContigVariant(const ContigVariant&) = delete; // Use deep copy.
  ~ContigVariant() = default;

  ContigVariant& operator=(const ContigVariant&) = delete; // Use deep copy.

  // Always use deep copy when modifying this object (filter and set operations).
  std::shared_ptr<ContigVariant> deepCopy() const;

  void addVariant(std::shared_ptr<const Variant>& variant_ptr);
  const ContigId_t& contigId() const { return contig_id_; }
  size_t variantCount() const { return offset_variant_map_.size(); }

  // Set functions.
  bool isElement(const Variant& variant) const;
  std::shared_ptr<ContigVariant> Union(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  std::shared_ptr<ContigVariant> Intersection(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  std::shared_ptr<ContigVariant> Difference(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;

  std::shared_ptr<ContigVariant> filterVariants(const VariantFilter& filter) const;

  const OffsetVariantMap& getMap() const { return offset_variant_map_; }

  size_t size() const { return offset_variant_map_.size(); }

  // All variants in [start, end) - note that end points past the last variant; end = (last + 1).
  bool getSortedVariants(ContigOffset_t start,
                         ContigOffset_t end,
                         OffsetVariantMap& variant_map) const;


private:

  ContigId_t contig_id_;
  OffsetVariantMap offset_variant_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  bool getCodingSortedVariants(ContigId_t contig_id,
                               ContigOffset_t start,
                               ContigOffset_t end,
                               OffsetVariantMap& variant_map) const;

  const Attributes& attributes() const { return attributes_; }
  Attributes& attributes() { return attributes_; }
  void attributes(const Attributes& attributes) { attributes_ = attributes; }

  bool mutantProtein( const ContigId_t& contig_id,
                      const FeatureIdent_t& gene_id,
                      const FeatureIdent_t& sequence_id,
                      const std::shared_ptr<const GenomeDatabase>& genome_db,
                      std::shared_ptr<AminoSequence>& amino_sequence) const;

  // The coding variants in the variant_map are used to mutate the dna_sequence.
  static bool mutateDNA(const OffsetVariantMap& variant_map,
                        const FeatureIdent_t& sequence_id,
                        std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr);

private:

  GenomeId_t genome_id_;
  GenomeVariantMap genome_variant_map_;
  Attributes attributes_;

};



}   // namespace genome
}   // namespace kellerberrin

std::ostream & operator<<(std::ostream &os, const kellerberrin::genome::GenomeVariant& genome_variant);


#endif //KGL_VARIANT_DB_H
