//
// Created by kellerberrin on 3/01/18.
//

#ifndef KGL_VARIANT_DB_CONTIG_H
#define KGL_VARIANT_DB_CONTIG_H


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

  bool addVariant(std::shared_ptr<const Variant>& variant_ptr);

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


}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_DB_CONTIG_H
