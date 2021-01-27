//
// Created by kellerberrin on 1/7/20.
//

#ifndef KGL_VARIANT_DB_UNPHASED_CONTIG_H
#define KGL_VARIANT_DB_UNPHASED_CONTIG_H


#include "kel_utility.h"
#include "kgl_variant.h"
#include "kgl_variant_filter.h"
#include "kgl_variant_db_offset.h"
#include "kgl_variant_mutation.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object holds variants for each contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetDBMap = std::map<ContigOffset_t, std::unique_ptr<OffsetDB>>;

class ContigDB {

public:

  explicit ContigDB(ContigId_t contig_id) : contig_id_(std::move(contig_id)) {}
  virtual ~ContigDB() = default;

  ContigDB(const ContigDB &) = delete;
  [[nodiscard]] ContigDB &operator=(const ContigDB &) = delete; // Use deep copy.

  // Use this to copy the object.
  [[nodiscard]] std::shared_ptr<ContigDB> deepCopy() const { return filterVariants(TrueFilter()); }

  [[nodiscard]] const ContigId_t &contigId() const { return contig_id_; }

  // Unconditionally adds a variant to the contig (unique or not).
  [[nodiscard]]  bool addVariant(const std::shared_ptr<const Variant> &variant_ptr);

  // Only adds the variant if it does not exist at the offset (only unique).
  // Phasing information is removed, this function cannot be specified with a haploid or diploid population.
  [[nodiscard]]  bool addUniqueUnphasedVariant(const std::shared_ptr<const Variant> &variant_ptr);

  [[nodiscard]]  size_t variantCount() const;

  [[nodiscard]] const OffsetDBMap &getMap() const { return contig_offset_map_; }

  // Create a filtered contig.
  [[nodiscard]] std::shared_ptr<ContigDB> filterVariants(const VariantFilter &filter) const;

  // Filter this contig (efficient for large databases).
  // Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
  std::pair<size_t, size_t> inSituFilter(const VariantFilter &filter);

  // Processes all variants in the contig with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)
  template<class Obj, typename Func> bool processAll(Obj& object, Func objFunc) const;
  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const ContigReference> &contig_db_ptr) const;

  // Optionally returns the variant pointer if found in this contig
  [[nodiscard]] std::optional<std::shared_ptr<const Variant>> findVariant(const Variant& variant) const;

  // Returns a contig containing all variants in this contig that match the template contig.
  [[nodiscard]] std::shared_ptr<ContigDB> findContig(const std::shared_ptr<const ContigDB>& template_contig) const;

  [[nodiscard]] std::optional<OffsetDBArray> findOffsetArray(ContigOffset_t offset) const;

  bool getSortedVariants( PhaseId_t phase,
                          ContigOffset_t start,
                          ContigOffset_t end,
                          OffsetVariantMap& variant_map) const;

  // Retrieves a contig subset in the offset range [begin, end)
  [[nodiscard]] std::shared_ptr<ContigDB> subset(ContigOffset_t start, ContigOffset_t end) const;

private:

  ContigId_t contig_id_;
  OffsetDBMap contig_offset_map_;

  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex add_variant_mutex_;

  void checkUpstreamDeletion(OffsetVariantMap& variant_map) const;

};


// General purpose genome processing template.
// Processes all variants in the contig with class Obj and Func = &(bool Obj::objFunc(const std::shared_ptr<const Variant>))
template<class Obj, typename Func>
bool ContigDB::processAll(Obj& object, Func objFunc)  const {

  for (auto const& [offset, offset_ptr] : getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      if (not (object.*objFunc)(variant_ptr)) {

        ExecEnv::log().error("ContigDB::processAll<Obj, Func>; Problem executing general purpose template function at offset: {}", offset);
        return false;

      }

    }

  }

  return true;

}






} // namespace



#endif //KGL_VARIANT_DB_UNPHASED_CONTIG_H
