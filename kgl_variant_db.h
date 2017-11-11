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
  ContigVariant(const ContigVariant&) = default;
  ~ContigVariant() = default;

  ContigVariant& operator=(const ContigVariant&) = default;

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
  GenomeVariant(const GenomeVariant&) = default;
  ~GenomeVariant() = default;

  GenomeVariant& operator=(const GenomeVariant& genome_variant) = default;

  const GenomeId_t& genomeId() const { return genome_id_; }
  void genomeId(const GenomeId_t& contig_id) { genome_id_ = contig_id; }

  size_t contigCount() const { return genome_variant_map_.size(); }

  bool addContigVariant(std::shared_ptr<ContigVariant>& contig_variant);
  bool getContigVariant(const ContigId_t& contig_id, std::shared_ptr<ContigVariant>& contig_variant);

  bool addVariant(std::shared_ptr<const Variant> variant);

  size_t size() const;

  std::shared_ptr<GenomeVariant> filterVariants(const VariantFilter& filter) const;

  const GenomeVariantMap& contigMap() const { return genome_variant_map_; }

  bool isElement(const Variant& variant) const;
  std::shared_ptr<GenomeVariant> Union(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Intersection(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Difference(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;

  static std::shared_ptr<GenomeVariant> emptyGenomeVariant(const GenomeId_t& genome_id,
                                                           const std::shared_ptr<const GenomeDatabase>& genome_db);

  std::shared_ptr<GenomeVariant>
  disaggregateCompoundVariants(const std::shared_ptr<const GenomeDatabase>& genome_db) const;

  std::string output(char field_delimiter, VariantOutputIndex output_index) const;
  bool outputCSV(const std::string& file_name, VariantOutputIndex output_index) const;

  const Attributes& attributes() const { return attributes_; }
  Attributes& attributes() { return attributes_; }

private:

  GenomeId_t genome_id_;
  GenomeVariantMap genome_variant_map_;
  Attributes attributes_;

};



}   // namespace genome
}   // namespace kellerberrin

std::ostream & operator<<(std::ostream &os, const kellerberrin::genome::GenomeVariant& genome_variant);


#endif //KGL_VARIANT_DB_H
