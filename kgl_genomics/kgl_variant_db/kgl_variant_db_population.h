   //
// Created by kellerberrin on 13/08/18.
//

#ifndef KGL_VARIANT_DB_POPULATION_H
#define KGL_VARIANT_DB_POPULATION_H


#include "kgl_variant_db_genome.h"
#include "kgl_data_file_type.h"

#include <map>
#include <mutex>


namespace kellerberrin::genome {   //  organization::project

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The population level contains genomes, which in turn contain contigs (chromosomes),
// which turn contains offsets, which in turn contain an array of variants.
// The DataDB contains further file information; phased or unphased, diploid or haploid etc, etc.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeDBMap = std::map<GenomeId_t, std::shared_ptr<GenomeDB>>;

// Used in the processAll_MT()) templates.
template<class ObjFunc>
using MemberGenomeFunc = bool (ObjFunc::*)(std::shared_ptr<const GenomeDB>, const std::shared_ptr<const Variant>&);
using GenomeProcessFunc = std::function<bool(std::shared_ptr<const GenomeDB>, const std::shared_ptr<const Variant>&)>;

class PopulationDB : public DataDB {

public:

  PopulationDB(PopulationId_t population_id, DataSourceEnum data_source) : DataDB(data_source), population_id_(std::move(population_id)) {}
  PopulationDB(const PopulationDB&) = delete; // Use deep copy.
  ~PopulationDB() override { clear(); }  // Experimental, may be quicker than relying on smart pointer reference counting.

  // Preferred to fileId().
  [[nodiscard]] const std::string& populationId() const { return population_id_; }
  [[nodiscard]] const std::string& fileId() const override { return populationId(); }

  PopulationDB& operator=(const PopulationDB&) = delete; // Use deep copy.

  // Use this to copy the object. Just the trivial 'TrueFilter()'.
  [[nodiscard]] std::shared_ptr<PopulationDB> deepCopy() const;

  // Use this to empty the object. Just the trivial 'FalseFilter()'.
  std::pair<size_t, size_t> clear();

  // Create the genome variant if it does not exist.
  [[nodiscard]] std::optional<std::shared_ptr<GenomeDB>> getCreateGenome(const GenomeId_t& genome_id);

  // Retrieve a genome
  [[nodiscard]] std::optional<std::shared_ptr<GenomeDB>> getGenome(const GenomeId_t& genome_id) const;

  // Total variants held in this population, not unique variants.
  [[nodiscard]] size_t variantCount() const;

  // Returns all the unique variants in the population using the variant HGVS signature to determine uniqueness.
  [[nodiscard]] std::map<std::string, std::shared_ptr<const Variant>> uniqueVariants() const;

  // Create an equivalent population that is canonical variants, SNP are represented by '1X', Deletes by '1MnD'
  // and Inserts by '1MnI'. The population structure is re-created and is not a shallow copy.
  [[nodiscard]] std::unique_ptr<PopulationDB> canonicalPopulation() const;

  // Returns a list of contigs in the population and a sum of all variants in each contig_ref_ptr.
  // Useful for processing contig_ref_ptr based VCF files (could be multi-threaded for speed).
  [[nodiscard]] std::map<ContigId_t , size_t> contigCount() const;

  // Creates a filtered copy of the population database.
  // This is multi-threaded across genomes to be time efficient for large databases.
  // Important, CPU and memory efficient, but returns a shallow copy of the population.
  // If the filtered view needs to be used in another program scope then use
  // deepCopy() or selfFilter() to create a permanent view.
  [[nodiscard]] std::unique_ptr<PopulationDB> viewFilter(const BaseFilter& filter) const;


  // Filters the actual (this) population database, multi-threaded to be efficient for large databases.
  // selfFilter returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that remain after filtering.
  std::pair<size_t, size_t> selfFilter(const BaseFilter& filter);

  // ReturnType the underlying genome map.
  [[nodiscard]] const GenomeDBMap& getMap() const { return genome_map_; }

  // Unconditionally adds a genome to the population, returns false if the genome already exists.
  bool addGenome(const std::shared_ptr<GenomeDB>& genome);

  // Deletes any empty Genomes, returns number deleted.
  size_t trimEmpty();
  // The opposite of the above. Ensures that all genomes have an identical number of contigs, even if empty.
  // Returns the number of contigs that each genome now has,
  size_t squareContigs();
  // Unconditionally add a variant to the population.
  // This function is thread safe for concurrent updates.
  // The population structure cannot be 'read' while it is being updated.
  // Note that if a genome id occurs twice (homozygous) in the genome vector it will be
  // added twice to the relevant genome indicating the presence of a homozygous variant.
  [[nodiscard]] bool addVariant( const std::shared_ptr<const Variant>& variant_ptr,
                                 const std::vector<GenomeId_t>& genome_vector);

  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the reference genome.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const GenomeReference>& genome_db) const;
  // Compress a population into a single genome. Done when generating aggregate variant statistics for a population.
  [[nodiscard]] std::shared_ptr<GenomeDB> compressPopulation() const;
  // Compress a population into a single genome of unique (only) variants. Removes any variant phasing information.
  // Used to convert a Diploid/Haploid population to an unphased single genome. Source populations are unchanged.
  [[nodiscard]] std::shared_ptr<PopulationDB> uniqueUnphasedGenome() const;
  // Get the Info header, get the field header object from the first variant in the population.
  // Careful, this implicitly assumes that all variants in the population have the same DataMemoryBlock this only true
  // Of populations generated by a single VCF file.
  [[nodiscard]] std::optional<std::shared_ptr<const InfoEvidenceHeader>> getVCFInfoEvidenceHeader() const;
  // Utility function, processes all variants in the population with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)->bool.
  template<class Obj> bool processAll(Obj& object, MemberVariantFunc<Obj> objFunc) const {

    VariantProcessFunc callable = std::bind_front(objFunc, &object);
    return processAll(callable);

  }
  // Implementation
  bool processAll(const VariantProcessFunc& procFunc) const;

  // Same as above but multi-threaded. Note that Func = &Obj::objFunc(std::shared_ptr<const GenomeDB>, const shared_ptr<const Variant>&)->bool.
  template<class Obj> bool processAll_MT(Obj& object, MemberGenomeFunc<Obj> objFunc) const {

    GenomeProcessFunc callable = std::bind_front(objFunc, &object);
    return processAll_MT(callable);

  }
  // Implementation
  bool processAll_MT(const GenomeProcessFunc& objFunc) const;

private:

  GenomeDBMap genome_map_;
  PopulationId_t population_id_;
  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex add_variant_mutex_;
  // mutex to lock the structure when performing a selfFilter.
  mutable std::mutex insitufilter_mutex_;

};


}   // end namespace




#endif //KGL_VARIANT_DB_POPULATION_H
