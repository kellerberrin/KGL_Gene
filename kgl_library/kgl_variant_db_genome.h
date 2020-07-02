//
// Created by kellerberrin on 23/04/18.
//

#ifndef KGL_VARIANT_DB_UNPHASED_H
#define KGL_VARIANT_DB_UNPHASED_H



#include "kel_utility.h"
#include "kgl_variant.h"
#include "kgl_variant_db_offset.h"
#include "kgl_variant_db_contig.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants indexed by genome.
// This object hold variants for each genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class VariantContig>
using ContigVariantMap = std::map<ContigId_t, std::shared_ptr<VariantContig>>;

template<class VariantContig>
class GenomeVariantArray;

using UnphasedGenome = GenomeVariantArray<UnphasedContig>;

template<class VariantContig>
class GenomeVariantArray {

public:

  explicit GenomeVariantArray(const GenomeId_t& genome_id) : genome_id_(genome_id) {}
  virtual ~GenomeVariantArray() = default;

  GenomeVariantArray(const GenomeVariantArray&) = delete;
  [[nodiscard]] GenomeVariantArray& operator=(const GenomeVariantArray&) = delete; // Use deep copy.

  [[nodiscard]] std::shared_ptr<GenomeVariantArray> deepCopy() const; // Use this to copy the object.

  // unconditionally merge (retains duplicates) genomes and variants into this genome.
  [[nodiscard]] size_t mergeGenome(const std::shared_ptr<const GenomeVariantArray>& merge_genome);

  // A copy of this genome that only contains unique variants (all duplicates removed).
  [[nodiscard]] std::shared_ptr<GenomeVariantArray> uniqueGenome() const;

  [[nodiscard]] size_t variantCount() const;

  [[nodiscard]] bool addVariant(const std::shared_ptr<const Variant>& variant);

  // The first bool is normal operation. The second bool is if a unique variant was added to the genome.
  [[nodiscard]] bool addUniqueVariant(const std::shared_ptr<const Variant>& variant);

  [[nodiscard]] const GenomeId_t& genomeId() const { return genome_id_; }

  [[nodiscard]] std::shared_ptr<GenomeVariantArray<VariantContig>> filterVariants(const VariantFilter& filter) const;

  [[nodiscard]] const ContigVariantMap<VariantContig>& getMap() const { return contig_map_; }

  [[nodiscard]] std::optional<std::shared_ptr<VariantContig>> getCreateContig(const ContigId_t& contig_id);
  // Processes all variants in the population with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)
  template<class Obj, typename Func> bool processAll(Obj& object, Func objFunc) const;
  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const GenomeReference>& genome_db_ptr) const;

private:

  ContigVariantMap<VariantContig> contig_map_;
  GenomeId_t genome_id_;

  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex add_variant_mutex_;

  [[nodiscard]] bool addContig(const std::shared_ptr<VariantContig>& contig_ptr);

};


// General purpose genome processing template.
// Processes all variants in the genome with class Obj and Func = &(bool Obj::objFunc(const std::shared_ptr<const Variant>))
template<class VariantContig>
template<class Obj, typename Func>
bool GenomeVariantArray<VariantContig>::processAll(Obj& object, Func objFunc)  const {

    for (auto const& contig : getMap()) {

      for (auto const& variant_vector : contig.second->getMap()) {

        for (auto const& variant_ptr : variant_vector.second->getVariantArray()) {

          if (not (object.*objFunc)(variant_ptr)) {

            ExecEnv::log().error("UnphasedPopulation::processAll<Obj, Func>; Problem executing general purpose template function.");
            return false;

          }

        }

      }

    }

  return true;

}


// Use this to copy the object.
template<class VariantContig>
std::shared_ptr<GenomeVariantArray<VariantContig>> GenomeVariantArray<VariantContig>::deepCopy() const {

  std::shared_ptr<GenomeVariantArray> genome_copy(std::make_shared<GenomeVariantArray>(genomeId()));

  for (auto const& [contig_id, contig_ptr] :  getMap()) {

    if (not genome_copy->addContig(contig_ptr->deepCopy())) {

      ExecEnv::log().critical("GenomeVariantArray::deepCopy(), Genome: {}, Unable to deepcopy Contig: {}", genomeId(), contig_id);

    }

  }

  return genome_copy;

}


// unconditionally merge (retains duplicates) genomes and variants into this genome.
template<class VariantContig>
size_t GenomeVariantArray<VariantContig>::mergeGenome(const std::shared_ptr<const GenomeVariantArray>& merge_genome) {

  if (not merge_genome->processAll(*this, &GenomeVariantArray::addVariant)) {

    ExecEnv::log().error("GenomeVariantArray::mergeGenome, problem merging genome: {} into genome: {}",
                         merge_genome->genomeId(), genomeId());

  }

  return variantCount();

}


template<class VariantContig>
std::shared_ptr<GenomeVariantArray<VariantContig>> GenomeVariantArray<VariantContig>::uniqueGenome() const {

  std::shared_ptr<GenomeVariantArray> unique_genome(std::make_shared<GenomeVariantArray>(genomeId()));

  if (not processAll(*unique_genome, &GenomeVariantArray::addUniqueVariant)) {

    ExecEnv::log().error("GenomeVariantArray::mergeGenome, problem creating unique variant genome from genome: {}", genomeId());

  }

  return unique_genome;

}


template<class VariantContig>
bool GenomeVariantArray<VariantContig>::addVariant(const std::shared_ptr<const Variant>& variant) {

  auto contig_opt = getCreateContig(variant->contigId());
  if (not contig_opt) {

    ExecEnv::log().error("GenomeVariantArray::addVariant(), Genome: {} could not get or create Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  if (not contig_opt.value()->addVariant(variant)) {

    ExecEnv::log().error("GenomeVariantArray::addVariant(), Genome: {} could not add variant to Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  return true;

}


// The first bool is normal operation. The second bool is if a unique variant was added to the genome.
template<class VariantContig>
bool GenomeVariantArray<VariantContig>::addUniqueVariant(const std::shared_ptr<const Variant>& variant) {

  auto contig_opt = getCreateContig(variant->contigId());
  if (not contig_opt) {

    ExecEnv::log().error("GenomeVariantArray::addVariant(), Genome: {} could not get or create Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  bool result = contig_opt.value()->addUniqueVariant(variant);
  if (not result) {

    ExecEnv::log().error("GenomeVariantArray::addVariant(), Genome: {} could not add variant to Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  return true;

}


template<class VariantContig>
std::optional<std::shared_ptr<VariantContig>> GenomeVariantArray<VariantContig>::getCreateContig(const ContigId_t& contig_id) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    return result->second;

  } else {

    auto contig_ptr = std::make_shared<VariantContig>(contig_id);
    auto insert_result = contig_map_.try_emplace(contig_id, contig_ptr);

    if (not insert_result.second) {

      ExecEnv::log().error("GenomeVariantArray::getCreateContig(), Could not add contig: {} to genome : {}", contig_id, genomeId());
      return std::nullopt;

    }

    return contig_ptr;

  }

}

template<class VariantContig>
bool GenomeVariantArray<VariantContig>::addContig(const std::shared_ptr<VariantContig>& contig_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_map_.try_emplace(contig_ptr->contigId(), contig_ptr);

  if (not result.second) {

    ExecEnv::log().error("GenomeVariantArray::addContig(); could not add contig: {} (duplicate) to the genome", contig_ptr->contigId());

  }

  return result.second;

}


template<class VariantContig>
size_t GenomeVariantArray<VariantContig>::variantCount() const {

  size_t variant_count = 0;

  for (auto contig : getMap()) {

    variant_count += contig.second->variantCount();

  }

  return variant_count;

}


template<class VariantContig>
std::shared_ptr<GenomeVariantArray<VariantContig>> GenomeVariantArray<VariantContig>::filterVariants(const VariantFilter& filter) const {

  std::shared_ptr<GenomeVariantArray> filtered_genome_ptr(std::make_shared<GenomeVariantArray>(genomeId()));

  for (const auto& [contig_id, contig_ptr] : getMap()) {

    auto filtered_contig = contig_ptr->filterVariants(filter);
    if (not filtered_genome_ptr->addContig(filtered_contig)) {

      ExecEnv::log().error("GenomeVariantArray::filterVariants(), Genome: {}, Unable to inserted filtered Contig: {}", genomeId(), contig_id);

    }

  }

  return filtered_genome_ptr;

}


// Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
// The second integer is the number variants that pass inspection by comparison to the genome database.
template<class VariantContig>
std::pair<size_t, size_t> GenomeVariantArray<VariantContig>::validate(const std::shared_ptr<const GenomeReference>& genome_db_ptr) const {

  std::pair<size_t, size_t> genome_count{0, 0};
  for (auto const& [contig_id, contig_ptr] : getMap()) {

    std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db_ptr->getContigSequence(contig_id);

    if (not contig_opt) {

      ExecEnv::log().error("GenomeVariantArray::validate, No matching contig found in GenomeDatabase for Variant Contig: {}", contig_id);
      continue;

    }

    std::pair<size_t, size_t> contig_count = contig_ptr->validate(contig_opt.value());

    if (contig_count.first != contig_count.second) {

      ExecEnv::log().warn("GenomeVariantArray::validate(), Genome: {} Validation Failed in Contig: {}, Total Variants: {} Validated: {}",
                          genomeId(), contig_id, contig_count.first, contig_count.second);

    }

    genome_count.first += contig_count.first;
    genome_count.second += contig_count.second;

  }

  return genome_count;

}


}   // end namespace


#endif //KGL_KGL_VARIANT_DB_UNPHASED_H
