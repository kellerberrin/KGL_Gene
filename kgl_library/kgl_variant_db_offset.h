//
// Created by kellerberrin on 1/7/20.
//

#ifndef KGL_VARIANT_DB_OFFSET_H
#define KGL_VARIANT_DB_OFFSET_H


#include "kgl_variant.h"

#include <vector>
#include <forward_list>
#include <memory>


namespace kellerberrin::genome {   //  organization level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the interface class for the low-level variant container.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using OffsetVariantArray = std::vector<std::shared_ptr<const Variant>>;

class UnphasedContigOffset {

public:

  UnphasedContigOffset() = default;
  virtual ~UnphasedContigOffset() = default;

  UnphasedContigOffset(const UnphasedContigOffset &) = delete;
  UnphasedContigOffset& operator=(const UnphasedContigOffset &) = delete;

  [[nodiscard]] virtual OffsetVariantArray getVariantArray() const = 0;

  virtual void addVariant(std::shared_ptr<const Variant> variant_ptr) = 0;

private:


};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants until they can be phased.
// This object holds an unphased vector of variants for each offset.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


class UnphasedContigVectorOffset : public UnphasedContigOffset {

public:

  UnphasedContigVectorOffset() = default;
  UnphasedContigVectorOffset(const UnphasedContigVectorOffset &) = delete;
  ~UnphasedContigVectorOffset() override = default;

  [[nodiscard]] OffsetVariantArray getVariantArray() const override { return variant_vector_; }

  void addVariant(std::shared_ptr<const Variant> variant_ptr) override { variant_vector_.emplace_back(variant_ptr); }

private:

  OffsetVariantArray variant_vector_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Variant collection implemented as a list.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


class UnphasedContigListOffset : public UnphasedContigOffset {

public:

  UnphasedContigListOffset() = default;
  UnphasedContigListOffset(const UnphasedContigListOffset &) = delete;
  ~UnphasedContigListOffset() override = default;

  [[nodiscard]] OffsetVariantArray getVariantArray() const override {

    OffsetVariantArray variant_vector;

    for (const auto& variant : variant_list_) {

      variant_vector.emplace_back(variant);

    }

    return variant_vector;

  }

  void addVariant(std::shared_ptr<const Variant> variant_ptr) override { variant_list_.emplace_front(variant_ptr); }

private:

  std::forward_list<std::shared_ptr<const Variant>> variant_list_;

};



////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Hold phased diploid variants (1000 Genome project).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


class DiploidOffset : public UnphasedContigOffset {

public:

  DiploidOffset() = default;
  DiploidOffset(const DiploidOffset &) = delete;
  ~DiploidOffset() override = default;

  [[nodiscard]] OffsetVariantArray getVariantArray() const override {

    OffsetVariantArray variant_vector;

    for (const auto& variant : variant_list_) {

      variant_vector.emplace_back(variant);

    }

    return variant_vector;

  }

  void addVariant(std::shared_ptr<const Variant> variant_ptr) override { variant_list_.emplace_front(variant_ptr); }

private:

  std::forward_list<std::shared_ptr<const Variant>> variant_list_;

};







} // namespace


#endif //KGL_KGL_VARIANT_DB_OFFSET_H
