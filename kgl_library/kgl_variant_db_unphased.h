//
// Created by kellerberrin on 23/04/18.
//

#ifndef KGL_VARIANT_DB_UNPHASED_H
#define KGL_VARIANT_DB_UNPHASED_H



#include "kel_utility.h"
#include "kgl_variant.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants until they can be phased.
// This object holds variants for each contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////



using UnphasedVectorVariantCount = std::vector<std::shared_ptr<const Variant>>;
using UnphasedOffsetMap = std::map<ContigOffset_t, UnphasedVectorVariantCount>;
class UnphasedContig {

public:

  explicit UnphasedContig(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  UnphasedContig(const UnphasedContig&) = delete;
  virtual ~UnphasedContig() = default;

  [[nodiscard]] UnphasedContig& operator=(const UnphasedContig&) = delete; // Use deep copy.

  [[nodiscard]] std::shared_ptr<UnphasedContig> deepCopy() const; // Use this to copy the object.

  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }

  // Unconditionally adds a variant to the contig (unique or not).
  [[nodiscard]]  bool addVariant(const std::shared_ptr<const Variant>& variant_ptr);

  // Test if an equivalent variant already exists in the contig.
  [[nodiscard]] bool variantExists(const std::shared_ptr<const Variant>& variant);
  // Only adds the variant if it does not already exist.
  [[nodiscard]] bool addUniqueVariant(const std::shared_ptr<const Variant>& variant);

  [[nodiscard]]  size_t variantCount() const;

  [[nodiscard]] const UnphasedOffsetMap& getMap() const { return contig_offset_map_; }


  [[nodiscard]] std::shared_ptr<UnphasedContig> filterVariants(const VariantFilter& filter) const;

  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const ContigReference>& contig_db_ptr) const;

private:

  ContigId_t contig_id_;
  UnphasedOffsetMap contig_offset_map_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants until they can be phased.
// This object hold variants for each genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using UnphasedContigMap = std::map<ContigId_t, std::shared_ptr<UnphasedContig>>;
class UnphasedGenome {

public:

  explicit UnphasedGenome(const GenomeId_t& genome_id) : genome_id_(genome_id) {}
  UnphasedGenome(const UnphasedGenome&) = delete;
  virtual ~UnphasedGenome() = default;

  [[nodiscard]] UnphasedGenome& operator=(const UnphasedGenome&) = delete; // Use deep copy.

  [[nodiscard]] std::shared_ptr<UnphasedGenome> deepCopy() const; // Use this to copy the object.

  // unconditionally merge (retains duplicates) genomes and variants into this genome.
  [[nodiscard]] size_t mergeGenome(const std::shared_ptr<const UnphasedGenome>& merge_genome);

  // A copy of this genome that only contains unique variants (all duplicates removed).
  [[nodiscard]] std::shared_ptr<UnphasedGenome> uniqueGenome() const;

  [[nodiscard]] size_t variantCount() const;

  [[nodiscard]] bool addVariant(const std::shared_ptr<const Variant>& variant);

  // The first bool is normal operation. The second bool is if a unique variant was added to the genome.
  [[nodiscard]] bool addUniqueVariant(const std::shared_ptr<const Variant>& variant);

  [[nodiscard]] const GenomeId_t& genomeId() const { return genome_id_; }

  [[nodiscard]] std::shared_ptr<UnphasedGenome> filterVariants(const VariantFilter& filter) const;

  [[nodiscard]] const UnphasedContigMap& getMap() const { return contig_map_; }

  [[nodiscard]] std::optional<std::shared_ptr<UnphasedContig>> getCreateContig(const ContigId_t& contig_id);
  // Processes all variants in the population with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)
  template<class Obj, typename Func> bool processAll(Obj& object, Func objFunc) const;
  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const GenomeReference>& genome_db_ptr) const;

private:

  UnphasedContigMap contig_map_;
  GenomeId_t genome_id_;

  [[nodiscard]] bool addContig(const std::shared_ptr<UnphasedContig>& contig_ptr);

};


// General purpose genome processing template.
// Processes all variants in the genome with class Obj and Func = &(bool Obj::objFunc(const std::shared_ptr<const Variant>))
template<class Obj, typename Func>
bool UnphasedGenome::processAll(Obj& object, Func objFunc)  const {

    for (auto const& contig : getMap()) {

      for (auto const& variant_vector : contig.second->getMap()) {

        for (auto const& variant_ptr : variant_vector.second) {

          if (not (object.*objFunc)(variant_ptr)) {

            ExecEnv::log().error("UnphasedPopulation::processAll<Obj, Func>; Problem executing general purpose template function.");
            return false;

          }

        }

      }

    }

  return true;

}




}   // end namespace


#endif //KGL_KGL_VARIANT_DB_UNPHASED_H
