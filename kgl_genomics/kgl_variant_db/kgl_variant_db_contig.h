//
// Created by kellerberrin on 1/7/20.
//

#ifndef KGL_VARIANT_DB_UNPHASED_CONTIG_H
#define KGL_VARIANT_DB_UNPHASED_CONTIG_H


#include "kel_utility.h"
#include "kgl_variant_db.h"
#include "kgl_variant_db_offset.h"
#include "kgl_mutation_variant.h"


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

  ContigDB(const ContigDB &) = delete; // Use deep copy.
  ContigDB &operator=(const ContigDB &) = delete; // Use deep copy.

  // Use this to copy the object.
  [[nodiscard]] std::shared_ptr<ContigDB> deepCopy() const;

  [[nodiscard]] const ContigId_t &contigId() const { return contig_id_; }

  // Unconditionally adds a variant to the contig (unique or not).
  [[nodiscard]]  bool addVariant(const std::shared_ptr<const Variant> &variant_ptr);


  [[nodiscard]]  size_t variantCount() const;

  [[nodiscard]] const OffsetDBMap &getMap() const { return contig_offset_map_; }

  // Return a filtered copy of the contig.
  // Important, returns a shallow copy of the contig - only use for CPU/memory efficiency.
  [[nodiscard]] std::unique_ptr<ContigDB> viewFilter(const BaseFilter &filter) const;
  // Filter this contig (efficient for large databases).
  // Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
  std::pair<size_t, size_t> selfFilter(const BaseFilter &filter);

  // Deletes any empty Offsets, returns number deleted.
  size_t trimEmpty();
  // Processes all variants in the contig with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)
  template<class Obj> bool processAll(Obj& object, MemberVariantFunc<Obj> objFunc) const;
  // Implementation
  bool processAll(const VariantProcessFunc& objFunc) const;

  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const ContigReference> &contig_db_ptr) const;

  // Optionally returns the variant pointer if found in this contig
  [[nodiscard]] std::optional<std::shared_ptr<const Variant>> findVariant(const Variant& variant) const;

  // Returns a contig containing all variants in this contig that match the template contig.
  [[nodiscard]] std::shared_ptr<ContigDB> findContig(const std::shared_ptr<const ContigDB>& template_contig) const;

  // Returns a variant offset array (if it exists) at a specified offset within the contig.
  [[nodiscard]] std::optional<OffsetDBArray> findOffsetArray(ContigOffset_t offset) const;

  bool getSortedVariants( VariantPhase phase,
                          ContigOffset_t start,
                          ContigOffset_t end,
                          OffsetVariantMap& variant_map) const;

  // Find the [lower, upper] offsets of a contig, {0,0} is returned if empty.
  [[nodiscard]] std::pair<ContigOffset_t, ContigOffset_t> offsetBounds() const;
  // Retrieves a contig subset in the offset range [begin, end)
  [[nodiscard]] std::shared_ptr<ContigDB> subset(ContigOffset_t start, ContigOffset_t end) const;

  // Unconditionally add all the variants in the supplied contig to this contig.
  bool merge(const std::shared_ptr<const ContigDB>& contig) { return contig->processAll(*this, &ContigDB::addVariant); }

  //Set Functions.
  // setIntersection returns a contig that contains variants present in both contigs.
  // The VariantEquality flag determines whether variant phase is used in the equality.
  [[nodiscard]] std::unique_ptr<ContigDB> setIntersection(const ContigDB& intersection_contig, VariantEquality variant_equality) const;
  // setComplement returns a contig that contains variants present in this contig but not present in the complement_contig.
  [[nodiscard]] std::unique_ptr<ContigDB> setComplement(const ContigDB& complement_contig, VariantEquality variant_equality) const;
  // setUnion returns a contig that contains the set union of variants present in this contig and the union_contig.
  // Specifically, identical variants, defined by the VariantEquality flag as including phase or not, are ignored.
  [[nodiscard]] std::unique_ptr<ContigDB> setUnion(const ContigDB& union_contig, VariantEquality variant_equality) const;

private:


  ContigId_t contig_id_;
  OffsetDBMap contig_offset_map_;

  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex lock_contig_mutex_;

  void checkUpstreamDeletion(OffsetVariantMap& variant_map) const;

  // Unconditionally adds an offset
  [[nodiscard]]  bool addOffset(ContigOffset_t offset, std::unique_ptr<OffsetDB> offset_db);

};


// General purpose genome processing template.
// Processes all variants in the contig with class Obj and Func = &(bool Obj::objFunc(const std::shared_ptr<const Variant>&))
template<class Obj>
bool ContigDB::processAll(Obj& object, MemberVariantFunc<Obj> objFunc)  const {

  for (auto const& [offset, offset_ptr] : getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      if (not std::invoke(objFunc, object, variant_ptr)) {

        ExecEnv::log().error("ContigDB::processAll<Obj, Func>; Problem executing general purpose template function at offset: {}", offset);
        return false;

      }

    }

  }

  return true;

}






} // namespace



#endif //KGL_VARIANT_DB_UNPHASED_CONTIG_H
