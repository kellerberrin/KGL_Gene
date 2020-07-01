//
// Created by kellerberrin on 1/7/20.
//

#ifndef KGL_VARIANT_DB_UNPHASED_CONTIG_H
#define KGL_VARIANT_DB_UNPHASED_CONTIG_H


#include "kel_utility.h"
#include "kgl_variant.h"
#include "kgl_variant_db_offset.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants until they can be phased.
// This object holds variants for each contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class VariantArray> class ContigOffsetVariant;

using UnphasedContig = ContigOffsetVariant<UnphasedContigListOffset>;

using UnphasedOffsetMap = std::map<ContigOffset_t, std::unique_ptr<UnphasedContigOffset>>;

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

  // Test if an equivalent variant already exists in the contig.
  [[nodiscard]] bool variantExists(const std::shared_ptr<const Variant> &variant);

  // Only adds the variant if it does not already exist.
  [[nodiscard]] bool addUniqueVariant(const std::shared_ptr<const Variant> &variant);

  [[nodiscard]]  size_t variantCount() const;

  [[nodiscard]] const UnphasedOffsetMap &getMap() const { return contig_offset_map_; }

  [[nodiscard]] std::shared_ptr<ContigOffsetVariant> filterVariants(const VariantFilter &filter) const;

  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the genome database.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const ContigReference> &contig_db_ptr) const;

private:

  ContigId_t contig_id_;
  UnphasedOffsetMap contig_offset_map_;

};


// Unconditionally adds a variant to the contig.
template<class VariantArray>
bool ContigOffsetVariant<VariantArray>::addVariant(const std::shared_ptr<const Variant> &variant_ptr) {

  auto result = contig_offset_map_.find(variant_ptr->offset());

  if (result != contig_offset_map_.end()) {
    // Variant offset exists.

    result->second->addVariant(variant_ptr);

  } else {
    // add the new offset.
    std::unique_ptr<UnphasedContigOffset> offset_array_ptr(std::make_unique<VariantArray>());
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

  std::shared_ptr<UnphasedContig> contig_copy(std::make_shared<UnphasedContig>(contigId()));

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


// Test if an equivalent variant already exists in the contig.
template<class VariantArray>
bool ContigOffsetVariant<VariantArray>::variantExists(const std::shared_ptr<const Variant> &variant) {

  auto result = contig_offset_map_.find(variant->offset());

  if (result != contig_offset_map_.end()) {
    // Variant offset exists.

    for (auto const &existingvariant : result->second->getVariantArray()) {

      if (existingvariant->equivalent(*variant)) {

        // found the variant.
        return true;

      }
    }
    // If we fall through the loop then no equivalent variant was found.
    return false;

  } else {

    // No variants found at the offset.
    return false;

  }

}

// Adds a unique variant to the contig. No addition if not unique.
template<class VariantArray>
bool ContigOffsetVariant<VariantArray>::addUniqueVariant(const std::shared_ptr<const Variant> &variant) {

  if (variantExists(variant)) {
// don't add.

    // Operation normal but variant not unique.
    return true;

  }

  return addVariant(variant);

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

  std::shared_ptr<UnphasedContig> filtered_contig_ptr(std::make_shared<UnphasedContig>(contigId()));

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

      if (contig_sequence_ptr->subSequence(variant_ptr->offset(), variant_ptr->reference().length()) ==
          variant_ptr->reference()) {

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


} // namespace



#endif //KGL_VARIANT_DB_UNPHASED_CONTIG_H
