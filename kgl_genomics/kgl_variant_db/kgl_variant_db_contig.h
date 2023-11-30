//
// Created by kellerberrin on 1/7/20.
//

#ifndef KGL_VARIANT_DB_CONTIG_H
#define KGL_VARIANT_DB_CONTIG_H


#include "kgl_variant_db_offset.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object holds variants for each contig_ref_ptr.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetDBMap = std::map<ContigOffset_t, std::shared_ptr<OffsetDB>>;

class ContigDB {

public:

  explicit ContigDB(ContigId_t contig_id) : contig_id_(std::move(contig_id)) {}
  virtual ~ContigDB() = default;

  ContigDB(const ContigDB &) = delete; // Use deep copy.
  ContigDB &operator=(const ContigDB &) = delete; // Use deep copy.

  // Use this to copy the object.
  [[nodiscard]] std::shared_ptr<ContigDB> deepCopy() const;

  [[nodiscard]] const ContigId_t &contigId() const { return contig_id_; }

  // Unconditionally adds a variant to the contig_ref_ptr (unique or not).
  [[nodiscard]]  bool addVariant(const std::shared_ptr<const Variant> &variant_ptr);


  [[nodiscard]]  size_t variantCount() const;

  [[nodiscard]] const OffsetDBMap &getMap() const { return contig_offset_map_; }

  // Return a filtered copy of the contig_ref_ptr.
  // Important, returns a shallow copy of the contig_ref_ptr - only use for CPU/memory efficiency.
  [[nodiscard]] std::unique_ptr<ContigDB> viewFilter(const BaseFilter &filter) const;
  // Filter this contig_ref_ptr (efficient for large databases).
  // Returns a std::pair with .first the reference number of variants, .second the filtered number of variants.
  std::pair<size_t, size_t> selfFilter(const BaseFilter &filter);

  // Deletes any empty Offsets, returns number deleted.
  size_t trimEmpty();
  // Processes all variants in the contig_ref_ptr with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)
  template<class Obj> bool processAll(Obj& object, MemberVariantFunc<Obj> objFunc) const;
  // Implementation
  bool processAll(const VariantProcessFunc& objFunc) const;

  // Create an equivalent contig_ref_ptr that has canonical variants, SNP are represented by '1X', Deletes by '1MnD'
  // and Inserts by '1MnI'. The contig_ref_ptr structure is re-created and is not a shallow copy.
  [[nodiscard]] std::unique_ptr<ContigDB> canonicalContig() const;

  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const ContigReference> &contig_db_ptr) const;

  // Returns a variant offset array (if it exists) at a specified offset within the contig_ref_ptr.
  [[nodiscard]] std::optional<OffsetDBArray> findOffsetArray(ContigOffset_t offset) const;

  // Unconditionally add all the variants in the supplied contig_ref_ptr to this contig.
  bool merge(const std::shared_ptr<const ContigDB>& contig) { return contig->processAll(*this, &ContigDB::addVariant); }

private:


  ContigId_t contig_id_;
  OffsetDBMap contig_offset_map_;

  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex lock_contig_mutex_;

  // Unconditionally adds an offset
  [[nodiscard]]  bool addOffset(ContigOffset_t offset, std::unique_ptr<OffsetDB> offset_db);

};


// General purpose genome processing template.
// Processes all variants in the contig_ref_ptr with class Obj and Func = &(bool Obj::objFunc(const std::shared_ptr<const Variant>&))
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



#endif //KGL_VARIANT_DB_CONTIG_H
