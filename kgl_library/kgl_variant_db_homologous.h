//
// Created by kellerberrin on 22/04/18.
//

#ifndef KGL_VARIANT_DB_HOMOLOGOUS_H
#define KGL_VARIANT_DB_HOMOLOGOUS_H



#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// HomologousVariant - All the variant features that map onto a homologous contig (haploid, diploid or polyploid)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetVariantMap = std::multimap<ContigOffset_t, std::shared_ptr<const Variant>>;
class HomologousVariant {

public:

  explicit HomologousVariant(const PhaseId_t& phase_id) : phase_id_(phase_id) {}
  HomologousVariant(const HomologousVariant&) = delete; // Use deep copy.
  ~HomologousVariant() = default;

  HomologousVariant& operator=(const HomologousVariant&) = delete; // Use deep copy.

  // Always use deep copy when modifying this object (filter and set operations).
  std::shared_ptr<HomologousVariant> deepCopy() const;

  bool addVariant(std::shared_ptr<const Variant>& variant_ptr);
  // Returns true if variant found and erased.
  bool eraseVariant(std::shared_ptr<const Variant>& variant_ptr);

  PhaseId_t phaseId() const { return phase_id_; }
  size_t variantCount() const { return offset_variant_map_.size(); }

  // Set functions.
  bool isElement(const Variant& variant) const;
  std::shared_ptr<HomologousVariant> Union(std::shared_ptr<const HomologousVariant> contig_variant_ptr) const;
  std::shared_ptr<HomologousVariant> Intersection(std::shared_ptr<const HomologousVariant> contig_variant_ptr) const;
  std::shared_ptr<HomologousVariant> Difference(std::shared_ptr<const HomologousVariant> contig_variant_ptr) const;

  std::shared_ptr<HomologousVariant> filterVariants(const VariantFilter& filter) const;

  const OffsetVariantMap& getMap() const { return offset_variant_map_; }

  size_t size() const { return offset_variant_map_.size(); }

  // All variants in [start, end) - note that end points past the last variant; end = (last + 1).
  bool getSortedVariants(ContigOffset_t start,
                         ContigOffset_t end,
                         OffsetVariantMap& variant_map) const;

private:

  PhaseId_t phase_id_;
  OffsetVariantMap offset_variant_map_;

};


}   // namespace genome
}   // namespace kellerberrin

// Not in kgl:: namespace.
std::ostream & operator<<(std::ostream &os, const kellerberrin::genome::OffsetVariantMap& variant_map);

#endif //KGL_KGL_VARIANT_DB_HOMOLOGOUS_H
