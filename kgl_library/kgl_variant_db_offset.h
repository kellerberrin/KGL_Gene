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

class VirtualContigOffset {

public:

  VirtualContigOffset() = default;
  virtual ~VirtualContigOffset() = default;

  VirtualContigOffset(const VirtualContigOffset &) = delete;
  VirtualContigOffset& operator=(const VirtualContigOffset &) = delete;

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


class UnphasedContigVectorOffset : public VirtualContigOffset {

public:

  UnphasedContigVectorOffset() = default;
  UnphasedContigVectorOffset(const UnphasedContigVectorOffset &) = delete;
  ~UnphasedContigVectorOffset() override = default;

  [[nodiscard]] OffsetVariantArray getVariantArray() const override {

    OffsetVariantArray variant_vector = variant_vector_;
    return variant_vector;

  }

  void addVariant(std::shared_ptr<const Variant> variant_ptr) override { variant_vector_.emplace_back(variant_ptr); }

private:

  OffsetVariantArray variant_vector_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Variant collection implemented as a list.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


class UnphasedContigListOffset : public VirtualContigOffset {

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
// Hold phased haploid variants (Pf3k).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

// There can be more than 1 variant per offset
// The VCF definition is that Indels are always specified with the reference
// as the nucleotide before the addition/deletion.
// Hence we can have an SNP and Indel occupy the same offset.
// However we issue a warning if more than 2 variants at an offset.

class HaploidOffset : public VirtualContigOffset {

public:

  HaploidOffset() {  variant_vector_.reserve(1); }
  HaploidOffset(const HaploidOffset &) = delete;
  ~HaploidOffset() override = default;

  [[nodiscard]] OffsetVariantArray getVariantArray() const override {

    OffsetVariantArray variant_vector = variant_vector_;

    if (variant_vector.empty()) {

      ExecEnv::log().error("DiploidOffset::getVariantArray, no variants found at location.");

    }

    return variant_vector;

  }

  void addVariant(std::shared_ptr<const Variant> variant_ptr) override {

    if (variant_ptr->phaseId() != Variant::HAPLOID_PHASED) {

      ExecEnv::log().error("HaploidOffset::addVariant, variant has incorrect phase: {}",
                           variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false));
    }

    variant_vector_.emplace_back(variant_ptr);

    if (variant_vector_.size() > 2) {

      ExecEnv::log().warn("HaploidOffset::addVariant, unexpected variant vector size: {}", variant_vector_.size());

      for (const auto& variant : variant_vector_) {

        ExecEnv::log().warn("HaploidOffset::addVariant, unexpected size variant: {}",
                             variant->output(' ', VariantOutputIndex::START_0_BASED, false));


      }

    }

  }

private:

  OffsetVariantArray variant_vector_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Hold phased diploid variants (1000 Genome project).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

// There can be more than 1 variant per offset
// The VCF definition is that Indels are always specified with the reference
// as the nucleotide before the addition/deletion.
// Hence we can have an SNP and Indel occupy the same offset.
// However we issue a warning if more than 2 variants at an offset.

class DiploidOffset : public VirtualContigOffset {

public:

  DiploidOffset() = default;
  DiploidOffset(const DiploidOffset &) = delete;
  ~DiploidOffset() override = default;

  [[nodiscard]] OffsetVariantArray getVariantArray() const override {

    OffsetVariantArray variant_vector;

    if (variant_A1_) {

      variant_vector.emplace_back(variant_A1_);

    }

    if (variant_A2_) {

      variant_vector.emplace_back(variant_A2_);

    }

    if (variant_B1_) {

      variant_vector.emplace_back(variant_B1_);

    }

    if (variant_B2_) {

      variant_vector.emplace_back(variant_B2_);

    }

    return variant_vector;

  }

  void addVariant(std::shared_ptr<const Variant> variant_ptr) override {

    if (variant_ptr->phaseId() == VariantSequence::DIPLOID_PHASE_A) {

      if (not variant_A1_) {

        variant_A1_ = variant_ptr;
        return;

      }

      if (not variant_A2_) {

        // Ignore identical variants.
        if (not variant_A1_->equivalent(*variant_ptr)) {

          variant_A2_ = variant_ptr;

        }

        return;

      }

      ExecEnv::log().warn("DiploidOffset::addVariant, unexpected A variant vector size: 3");
      ExecEnv::log().warn("DiploidOffset::addVariant, unexpected size variant: {}, VCF record: {}",
                          variant_A1_->output(' ', VariantOutputIndex::START_0_BASED, false),
                          variant_ptr->evidence().vcfRecordCount());
      ExecEnv::log().warn("DiploidOffset::addVariant, unexpected size variant: {}, VCF record: {}",
                          variant_A2_->output(' ', VariantOutputIndex::START_0_BASED, false),
                          variant_ptr->evidence().vcfRecordCount());
      ExecEnv::log().warn("DiploidOffset::addVariant, unexpected size variant: {}, VCF record: {}",
                          variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false),
                          variant_ptr->evidence().vcfRecordCount());

    }

    if (variant_ptr->phaseId() == VariantSequence::DIPLOID_PHASE_B) {

      if (not variant_B1_) {

        variant_B1_ = variant_ptr;
        return;

      }

      if (not variant_B2_) {

        // Ignore identical variants.
        if (not variant_B1_->equivalent(*variant_ptr)) {

          variant_B2_ = variant_ptr;

        }

        return;

      }

      ExecEnv::log().warn("DiploidOffset::addVariant, unexpected B variant vector size: 3");
      ExecEnv::log().warn("DiploidOffset::addVariant, unexpected size variant: {}, VCF record: {}",
                          variant_B1_->output(' ', VariantOutputIndex::START_0_BASED, false),
                          variant_ptr->evidence().vcfRecordCount());
      ExecEnv::log().warn("DiploidOffset::addVariant, unexpected size variant: {}, VCF record: {}",
                          variant_B2_->output(' ', VariantOutputIndex::START_0_BASED, false),
                          variant_ptr->evidence().vcfRecordCount());
      ExecEnv::log().warn("DiploidOffset::addVariant, unexpected size variant: {}, VCF record: {}",
                          variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false),
                          variant_ptr->evidence().vcfRecordCount());

    }

    if (variant_ptr->phaseId() != VariantSequence::DIPLOID_PHASE_A
        and variant_ptr->phaseId() != VariantSequence::DIPLOID_PHASE_B) {

      ExecEnv::log().error("DiploidOffset::addVariant, variant is not diploid phased: {} VCF record: {}",
                           variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false),
                           variant_ptr->evidence().vcfRecordCount());

    } else {

      ExecEnv::log().warn("DiploidOffset::addVariant, variant not processed: {}, VCF record: {}",
                          variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false),
                          variant_ptr->evidence().vcfRecordCount());

    }

  }

private:

  // There can be more than 1 variant per offset
  // The VCF definition is that Indels are always specified with the reference
  // as the nucleotide before the addition/deletion.
  // Hence we can have an SNP and Indel occupy the same offset.
  // However we issue a warning if more than 2 variants at an offset.

  std::shared_ptr<const Variant> variant_A1_;
  std::shared_ptr<const Variant> variant_A2_;
  std::shared_ptr<const Variant> variant_B1_;
  std::shared_ptr<const Variant> variant_B2_;

};







} // namespace


#endif //KGL_KGL_VARIANT_DB_OFFSET_H
