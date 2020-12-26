//
// Created by kellerberrin on 1/7/20.
//

#ifndef KGL_VARIANT_DB_UNPHASED_CONTIG_H
#define KGL_VARIANT_DB_UNPHASED_CONTIG_H


#include "kel_utility.h"
#include "kgl_variant.h"
#include "kgl_variant_db_offset.h"
#include "kgl_variant_mutation.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants until they can be phased.
// This object holds variants for each contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


using UnphasedOffsetMap = std::map<ContigOffset_t, std::unique_ptr<VariantArray>>;

class ContigOffsetVariant {

public:

  explicit ContigOffsetVariant(ContigId_t contig_id) : contig_id_(std::move(contig_id)) {}
  virtual ~ContigOffsetVariant() = default;

  ContigOffsetVariant(const ContigOffsetVariant &) = delete;
  [[nodiscard]] ContigOffsetVariant &operator=(const ContigOffsetVariant &) = delete; // Use deep copy.

  [[nodiscard]] std::shared_ptr<ContigOffsetVariant> deepCopy() const; // Use this to copy the object.

  [[nodiscard]] const ContigId_t &contigId() const { return contig_id_; }

  // Unconditionally adds a variant to the contig (unique or not).
  [[nodiscard]]  bool addVariant(const std::shared_ptr<const Variant> &variant_ptr);

  // Only adds the variant if it does not exist at the offset (only unique).
  // Phasing information is removed, this function cannot be specified with a haploid or diploid population.
  [[nodiscard]]  bool addUniqueUnphasedVariant(const std::shared_ptr<const Variant> &variant_ptr);

  [[nodiscard]]  size_t variantCount() const;

  [[nodiscard]] const UnphasedOffsetMap &getMap() const { return contig_offset_map_; }

  // Create a filtered contig.
  [[nodiscard]] std::shared_ptr<ContigOffsetVariant> filterVariants(const VariantFilter &filter) const;

  // Filter this contig (efficient for large databases).
  // Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
  std::pair<size_t, size_t> inSituFilter(const VariantFilter &filter);

  // Processes all variants in the contig with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)
  template<class Obj, typename Func> bool processAll(Obj& object, Func objFunc) const;
  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const ContigReference> &contig_db_ptr) const;

  [[nodiscard]] std::optional<std::shared_ptr<const Variant>> findVariant(const Variant& variant) const;

  [[nodiscard]] std::optional<OffsetVariantArray> findOffsetArray(ContigOffset_t offset) const;

  [[nodiscard]] bool getSortedVariants( PhaseId_t phase,
                                        ContigOffset_t start,
                                        ContigOffset_t end,
                                        OffsetVariantMap& variant_map) const;

private:

  ContigId_t contig_id_;
  UnphasedOffsetMap contig_offset_map_;

  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex add_variant_mutex_;

  void checkUpstreamDeletion(OffsetVariantMap& variant_map) const;

};


// General purpose genome processing template.
// Processes all variants in the contig with class Obj and Func = &(bool Obj::objFunc(const std::shared_ptr<const Variant>))
template<class Obj, typename Func>
bool ContigOffsetVariant::processAll(Obj& object, Func objFunc)  const {

  for (auto const& [offset, offset_ptr] : getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      if (not (object.*objFunc)(variant_ptr)) {

        ExecEnv::log().error("ContigOffsetVariant::processAll<Obj, Func>; Problem executing general purpose template function at offset: {}", offset);
        return false;

      }

    }

  }

  return true;

}






} // namespace



#endif //KGL_VARIANT_DB_UNPHASED_CONTIG_H
