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
#include "kgl_variant_db_homologous.h"
#include "kgl_variant.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using PloidyVector = std::vector<std::shared_ptr<HomologousVariant>>;
class ContigVariant {

public:

  explicit ContigVariant(const ContigId_t& contig_id, PhaseId_t phases)
  : contig_id_(contig_id) {

    for (PhaseId_t idx = 0; idx < phases; idx++) {

      std::shared_ptr<HomologousVariant> homologous_ptr(std::make_shared<HomologousVariant>(idx));
      ploidy_vector_.push_back(homologous_ptr);

    }

  }
  ContigVariant(const ContigVariant&) = delete; // Use deep copy.
  ~ContigVariant() = default;

  [[nodiscard]] ContigVariant& operator=(const ContigVariant&) = delete; // Use deep copy.

  // Always use deep copy when modifying this object (filter and set operations).
  [[nodiscard]] std::shared_ptr<ContigVariant> deepCopy() const;

  [[nodiscard]] bool addVariant(std::shared_ptr<const Variant>& variant_ptr);
  // Returns true if variant found and erased.
  [[nodiscard]] bool eraseVariant(std::shared_ptr<const Variant>& variant_ptr);

  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  [[nodiscard]] PhaseId_t ploidy() const { return static_cast<PhaseId_t>(ploidy_vector_.size()); }

  [[nodiscard]] size_t variantCount() const;

  // Set functions.
  [[nodiscard]] bool isElement(const Variant& variant) const;
  [[nodiscard]] std::shared_ptr<ContigVariant> Union(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  [[nodiscard]] std::shared_ptr<ContigVariant> Intersection(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  [[nodiscard]] std::shared_ptr<ContigVariant> Difference(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;

  [[nodiscard]] std::shared_ptr<ContigVariant> filterVariants(const VariantFilter& filter) const;

  [[nodiscard]] const PloidyVector& getVector() const { return ploidy_vector_; }

  [[nodiscard]] std::shared_ptr<const HomologousVariant> getHomologous(size_t index) const { return ploidy_vector_[index]; }

  // All variants in [start, end) - note that end points past the last variant; end = (last + 1).
  [[nodiscard]] bool getSortedVariants( PhaseId_t phase,
                                        ContigOffset_t start,
                                        ContigOffset_t end,
                                        OffsetVariantMap& variant_map) const;

  static constexpr PhaseId_t HAPLOID_HOMOLOGOUS_INDEX = 0;  // The first (and only) haploid homologous contig.

private:

  ContigId_t contig_id_;
  PloidyVector ploidy_vector_;

};


}   // end namespace

// Not in namespace.
[[nodiscard]] std::ostream & operator<<(std::ostream &os, const kellerberrin::genome::OffsetVariantMap& variant_map);


#endif //KGL_VARIANT_DB_CONTIG_H
