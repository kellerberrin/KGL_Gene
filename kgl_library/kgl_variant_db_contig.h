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


using UnphasedOffsetMap = std::map<ContigOffset_t, std::unique_ptr<VirtualContigOffset>>;

template<class VariantArray>
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

  [[nodiscard]]  size_t variantCount() const;

  [[nodiscard]] const UnphasedOffsetMap &getMap() const { return contig_offset_map_; }

  [[nodiscard]] std::shared_ptr<ContigOffsetVariant> filterVariants(const VariantFilter &filter) const;

  // Processes all variants in the contig with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)
  template<class Obj, typename Func> bool processAll(Obj& object, Func objFunc) const;
  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const ContigReference> &contig_db_ptr) const;

  [[nodiscard]] std::optional<std::shared_ptr<const Variant>> findVariant(const Variant& variant);

  [[nodiscard]] bool getSortedVariants( PhaseId_t phase,
                                        ContigOffset_t start,
                                        ContigOffset_t end,
                                        OffsetVariantMap& variant_map) const;

  constexpr static const PhaseId_t HAPLOID_HOMOLOGOUS_INDEX = VariantSequence::HAPLOID_PHASED;

private:

  ContigId_t contig_id_;
  UnphasedOffsetMap contig_offset_map_;

  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex add_variant_mutex_;

  void checkUpstreamDeletion(OffsetVariantMap& variant_map) const;

};


// General purpose genome processing template.
// Processes all variants in the contig with class Obj and Func = &(bool Obj::objFunc(const std::shared_ptr<const Variant>))
template<class VariantArray>
template<class Obj, typename Func>
bool ContigOffsetVariant<VariantArray>::processAll(Obj& object, Func objFunc)  const {

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


// Unconditionally adds a variant to the contig.
template<class VariantArray>
bool ContigOffsetVariant<VariantArray>::addVariant(const std::shared_ptr<const Variant> &variant_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_offset_map_.find(variant_ptr->offset());

  if (result != contig_offset_map_.end()) {
    // Variant offset exists.

    result->second->addVariant(variant_ptr);

  } else {
    // add the new offset.
    std::unique_ptr<VirtualContigOffset> offset_array_ptr(std::make_unique<VariantArray>());
    offset_array_ptr->addVariant(variant_ptr);
    auto insert_result = contig_offset_map_.try_emplace(variant_ptr->offset(), std::move(offset_array_ptr));

    if (not insert_result.second) {

      ExecEnv::log().error("UnphasedContig::addVariant(); Could not add variant offset: {} to the genome",
                           variant_ptr->offset());
      return false;

    }

  }

  return true;

}


// Use this to copy the object.
template<class VariantArray>
std::shared_ptr<ContigOffsetVariant<VariantArray>> ContigOffsetVariant<VariantArray>::deepCopy() const {

  std::shared_ptr<ContigOffsetVariant<VariantArray>> contig_copy(std::make_shared<ContigOffsetVariant<VariantArray>>(contigId()));

  for (auto const&[offset, variant_array] : getMap()) {

    for (auto const &variant_count : variant_array->getVariantArray()) {

      if (not contig_copy->addVariant(variant_count)) {

        ExecEnv::log().error("UnphasedContig::deepCopy; Cannot add Variant to Contig Copy : {}, at Offset: {}",
                             contig_copy->contigId(), offset);

      }

    }

  }

  return contig_copy;

}



// Counts the variants in a contug.
template<class VariantArray>
size_t ContigOffsetVariant<VariantArray>::variantCount() const {


  size_t variant_count = 0;

  for (auto const &offset_variant_vector : getMap()) {

    variant_count += offset_variant_vector.second->getVariantArray().size();

  }

  return variant_count;

}


// Creates a copy of the contig that only contains variants passing the filter condition.
template<class VariantArray>
std::shared_ptr<ContigOffsetVariant<VariantArray>> ContigOffsetVariant<VariantArray>::filterVariants(const VariantFilter &filter) const {

  std::shared_ptr<ContigOffsetVariant<VariantArray>> filtered_contig_ptr(std::make_shared<ContigOffsetVariant<VariantArray>>(contigId()));

  for (auto const&[offset, variant_vector] : getMap()) {

    for (auto const &variant_ptr : variant_vector->getVariantArray()) {

      if (filter.applyFilter(*variant_ptr)) {

        if (not filtered_contig_ptr->addVariant(variant_ptr)) {

          ExecEnv::log().error("UnphasedContig::filterVariants; Problem adding variant at offset: {}, to contig: {}",
                               offset, contigId());

        }

      }

    }

  }

  return filtered_contig_ptr;

}


template<class VariantArray>
[[nodiscard]] std::optional<std::shared_ptr<const Variant>> ContigOffsetVariant<VariantArray>::findVariant(const Variant& variant) {

  auto result = contig_offset_map_.find(variant.offset());

  if (result != contig_offset_map_.end()) {

    OffsetVariantArray variant_array = result->second->getVariantArray();

    for (auto const& variant_ptr : variant_array) {

      if (variant.homozygous(*variant_ptr)) {

        return variant_ptr;

      }

    }

    return std::nullopt;


  } else {

    return std::nullopt;

  }

}




template<class VariantArray>
std::pair<size_t, size_t> ContigOffsetVariant<VariantArray>::validate(const std::shared_ptr<const ContigReference> &contig_db_ptr) const {

  std::pair<size_t, size_t> contig_count{0, 0};

  std::shared_ptr<const DNA5SequenceContig> contig_sequence_ptr = contig_db_ptr->sequence_ptr();

  for (auto const&[offset, variant_vector] : getMap()) {

    contig_count.first += variant_vector->getVariantArray().size();

    if (offset >= contig_sequence_ptr->length()) {

      ExecEnv::log().error("UnphasedContig::validate,  Variant offset: {} exceeds total contig: {} size: {}", offset,
                           contig_db_ptr->contigId(), contig_sequence_ptr->length());
      continue;

    }

    for (auto const &variant_ptr : variant_vector->getVariantArray()) {

      if (not variant_ptr) {

        ExecEnv::log().error("UnphasedContig::validate, Unknown variant: {}",
                             variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false));
        continue;

      }

      if (contig_sequence_ptr->subSequence(variant_ptr->offset(), variant_ptr->reference().length()) == variant_ptr->reference()) {

        ++contig_count.second;

      } else {

        ExecEnv::log().error(
        "UnphasedContig::validate, Mismatch, at Contig Offset: {} Sequence is: {}, Variant Reference Sequence is: {}",
        variant_ptr->offset(),
        contig_sequence_ptr->subSequence(variant_ptr->offset(),
                                         variant_ptr->reference().length()).getSequenceAsString(),
        variant_ptr->reference().getSequenceAsString());

      }

    }

  }

  return contig_count;

}

template<class VariantArray>
bool ContigOffsetVariant<VariantArray>::getSortedVariants( PhaseId_t phase,
                                                           ContigOffset_t start,
                                                           ContigOffset_t end,
                                                           OffsetVariantMap& variant_map) const {


  auto lower_bound = contig_offset_map_.lower_bound(start);
  auto upper_bound = contig_offset_map_.upper_bound(end-1); //  [start, end)

  // If there is a prior variant that overlaps the start address, then push this onto the variant map.

  if (lower_bound != contig_offset_map_.end() and lower_bound != contig_offset_map_.begin()) {

    auto previous_offset_ptr = std::prev(lower_bound);

    OffsetVariantArray previous_offset_variants = previous_offset_ptr->second->getVariantArray();

    for (const auto& variant : previous_offset_variants) {

      if (variant->phaseId() == phase) {

        if (variant->offset() + variant->referenceSize() > start) {

          variant_map.emplace(variant->offset(), variant);

        }

      }

    }

  }



  for (auto it = lower_bound; it != upper_bound; ++it) {


    OffsetVariantArray previous_offset_variants = it->second->getVariantArray();

    for (const auto& variant : previous_offset_variants) {

      if (variant->phaseId() == phase) {

        variant_map.emplace(it->first, variant);

      }

    }

  }

  checkUpstreamDeletion(variant_map);

  return true;

}

template<class VariantArray>
void ContigOffsetVariant<VariantArray>::checkUpstreamDeletion(OffsetVariantMap& variant_map) const {

  for (auto iter = variant_map.begin(); iter != variant_map.end(); ++iter) {

    if (iter == variant_map.begin()) continue;

    int delete_size = std::prev(iter)->second->referenceSize() - std::prev(iter)->second->alternateSize();

    int offset_gap = iter->second->offset() - std::prev(iter)->second->offset();

    if (delete_size >= offset_gap) {

      ExecEnv::log().vinfo("checkUpstreamDeletion(), Upstream deletion detected: {}",
                           std::prev(iter)->second->output(' ', VariantOutputIndex::START_0_BASED, true));

      ExecEnv::log().vinfo("checkUpstreamDeletion(), Downstream variant removed from mutation: {}",
                           iter->second->output(' ', VariantOutputIndex::START_0_BASED, true));

      iter = variant_map.erase(iter);

      iter = std::prev(iter); // reset back to the previous iterator.

    }

  }

}






} // namespace



#endif //KGL_VARIANT_DB_UNPHASED_CONTIG_H
