//
// Created by kellerberrin on 23/04/18.
//

#ifndef KGL_VARIANT_DB_GENOME_H
#define KGL_VARIANT_DB_GENOME_H


#include "kgl_variant_db_contig.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object holds variants for each genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


using ContigDBMap = std::map<ContigId_t, std::shared_ptr<ContigDB>>;

class GenomeDB {

public:

  explicit GenomeDB(GenomeId_t genome_id) : genome_id_(std::move(genome_id)) {}
  virtual ~GenomeDB() = default;

  GenomeDB(const GenomeDB&) = delete; // Use deep copy.
  GenomeDB& operator=(const GenomeDB&) = delete; // Use deep copy.

  // Use this to copy the object.
  [[nodiscard]] std::shared_ptr<GenomeDB> deepCopy() const;

  // Unconditionally merge (retains duplicates) genomes and variants into this genome.
  [[nodiscard]] size_t mergeGenome(const std::shared_ptr<const GenomeDB>& merge_genome);

  [[nodiscard]] size_t variantCount() const;

  [[nodiscard]] bool addVariant(const std::shared_ptr<const Variant>& variant);

  [[nodiscard]] const GenomeId_t& genomeId() const { return genome_id_; }

  // Return a filtered copy of the genome.
  // Important, returns a shallow copy of the genome - only use for CPU/memory efficiency.
  [[nodiscard]] std::unique_ptr<GenomeDB> viewFilter(const BaseFilter& filter) const;
  // Filter this genome in Situ. (efficient for large databases).
  // Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
  std::pair<size_t, size_t> selfFilter(const BaseFilter &filter);

  // Deletes any empty Contigs, returns number deleted.
  size_t trimEmpty();

  // Create an equivalent genome that has canonical variants, SNP are represented by '1X', Deletes by '1MnD'
  // and Inserts by '1MnI'. The genome structure is re-created and is not a shallow copy.
  [[nodiscard]] std::unique_ptr<GenomeDB> canonicalGenome() const;

  // Get the underlying contig map.
  [[nodiscard]] const ContigDBMap& getMap() const { return contig_map_; }

  // Creates the contig if it does not exist.
  [[nodiscard]] std::optional<std::shared_ptr<ContigDB>> getCreateContig(const ContigId_t& contig_id);
  // Const version, Returns nullopt if the contig does not exist.
  [[nodiscard]] std::optional<std::shared_ptr<const ContigDB>> getContig(const ContigId_t& contig_id) const;
  // Returns nullopt if the contig does not exist.
  [[nodiscard]] std::optional<std::shared_ptr<ContigDB>> getContig(const ContigId_t& contig_id);
  // Processes all variants in the genome with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)
  template<class Obj> bool processAll(Obj& object, MemberVariantFunc<Obj> objFunc) const {

    VariantProcessFunc callable = std::bind_front(objFunc, &object);
    return processAll(callable);

  }
  // Implementation
  bool processAll(const VariantProcessFunc& procFunc) const;

  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const GenomeReference>& genome_db_ptr) const;

private:

  ContigDBMap contig_map_;
  GenomeId_t genome_id_;

  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex add_variant_mutex_;

  [[nodiscard]] bool addContig(std::shared_ptr<ContigDB> contig_ptr);

};



}   // end namespace


#endif //KGL_KGL_VARIANT_DB_UNPHASED_H
