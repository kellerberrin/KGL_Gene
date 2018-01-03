//
// Created by kellerberrin on 3/01/18.
//


#include <memory>
#include "kgl_variant_single.h"
#include "kgl_variant_compound.h"
#include "kgl_sequence_virtual_compare.h"
#include "kgl_variant_mutation.h"

namespace kgl = kellerberrin::genome;


// Mutate the DNA sequence using SNP variants
bool kgl::VariantMutation::mutateSNPs(const OffsetVariantMap& snp_variant_map,
                                      ContigOffset_t contig_offset,
                                      std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr) {

  // Mutate the base sequence.
  for (const auto& variant : snp_variant_map) {

    if (variant.second->isSingle()) {

      if (not mutateSingleSNP(variant.second, contig_offset, dna_sequence_ptr)) {

        ExecEnv::log().error("mutateSNPs(), problem mutating sequence with single snp: {}",
                             variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
        return false;

      }

    } else { // compound SNP.

      std::shared_ptr<const CompoundSNP> cmp_snp_ptr = std::dynamic_pointer_cast<const CompoundSNP>(variant.second);

      if (not cmp_snp_ptr) {

        ExecEnv::log().error("mutateSNPs(), should be CompoundSNP; unexpected variant: {}",
                             variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      } else {

        for (auto variant : cmp_snp_ptr->getMap()) {

          if (not mutateSingleSNP(variant.second, contig_offset, dna_sequence_ptr)) {

            ExecEnv::log().error("mutateSNPs(), problem mutating sequence with compound snp: {}",
                                 cmp_snp_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
            return false;
          }

        }

      }

    }

  }

  return true;

}

// Mutate the DNA sequence using SNP variants
bool kgl::VariantMutation::mutateSingleSNP(std::shared_ptr<const Variant> variant_ptr,
                                           ContigOffset_t contig_offset,
                                           std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr) {

  std::shared_ptr<const SNPVariant> snp_ptr = std::dynamic_pointer_cast<const SNPVariant>(variant_ptr);

  if (not snp_ptr) {

    ExecEnv::log().error("mutateSingleSNP, should be SNPVariant; unexpected variant: {}",
                         variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;

  }

  // Adjust the offset
  SignedOffset_t adjusted_offset = snp_ptr->offset() - contig_offset;

  // Check the offset
  if (adjusted_offset < 0 or adjusted_offset >= dna_sequence_ptr->length()) {

    ExecEnv::log().error("mutateSingleSNP(), calculated sequence offset: {} is out of range for sequence size: {}",
                         adjusted_offset, dna_sequence_ptr->length(),
                         snp_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;
  }

  ContigOffset_t sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);

  // Check the reference.
  if (snp_ptr->reference() != dna_sequence_ptr->at(sequence_offset)) {

    ExecEnv::log().warn("mutateSingleSNP(), reference base: {} does not match sequence base: {} at sequence offset: {}",
                        DNA5::convertToChar(snp_ptr->reference()),
                        DNA5::convertToChar(dna_sequence_ptr->at(sequence_offset)),
                        sequence_offset);

  }

  // Mutate the sequence
  dna_sequence_ptr->modifyBase(sequence_offset, snp_ptr->mutant());

  return true;

}

// Mutate the DNA sequence using Delete variants
bool kgl::VariantMutation::mutateDeletes(const OffsetVariantMap& snp_variant_map,
                                         ContigOffset_t contig_offset,
                                         std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                                         IndelAccountingMap& indel_accounting_map) {

  // Debug ***********
  return true;

  // Process the deletes from highest contig offset to lowest.
  // This avoids the need to adjust the contig offset of each delete.
  for (auto rit = snp_variant_map.rbegin(); rit !=  snp_variant_map.rend(); ++rit) {

    if (rit->second->isSingle()) {

      if (not mutateSingleDelete(rit->second, contig_offset, dna_sequence_ptr, indel_accounting_map)) {

        ExecEnv::log().error("mutateDeletes(), problem mutating sequence with single delete: {}",
                             rit->second->output(' ', VariantOutputIndex::START_0_BASED, true));
        return false;

      }

    } else { // compound SNP.

      std::shared_ptr<const CompoundDelete> cmp_del_ptr = std::dynamic_pointer_cast<const CompoundDelete>(rit->second);

      if (not cmp_del_ptr) {

        ExecEnv::log().error("mutateDeletes(), should be CompoundDelete; unexpected variant: {}",
                             rit->second->output(' ', VariantOutputIndex::START_0_BASED, true));

      } else {

        for (auto variant : cmp_del_ptr->getMap()) {

          if (not mutateSingleSNP(variant.second, contig_offset, dna_sequence_ptr)) {

            ExecEnv::log().error("mutateDeletes(), problem mutating sequence with compound snp: {}",
                                 cmp_del_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
            return false;
          }

        }

      }

    }

  }

  return true;

}


bool kgl::VariantMutation::mutateSingleDelete(std::shared_ptr<const Variant> variant_ptr,
                                              ContigOffset_t contig_offset,
                                              std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                                              IndelAccountingMap& indel_accounting_map) {

  std::shared_ptr<const DeleteVariant> del_ptr = std::dynamic_pointer_cast<const DeleteVariant>(variant_ptr);

  if (not del_ptr) {

    ExecEnv::log().error("mutateSingleDelete(), should be DeleteVariant; unexpected variant: {}",
                         variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;

  }


  // Adjust the offset
  SignedOffset_t adjusted_offset = del_ptr->offset() - contig_offset;

  adjusted_offset += adjustIndelOffsets(variant_ptr->offset(), indel_accounting_map);

  // Check the offset
  if (adjusted_offset < 0 or adjusted_offset >= dna_sequence_ptr->length()) {

    ExecEnv::log().error("mutateSingleDelete(), calculated sequence offset: {} is out of range for sequence size: {}",
                         adjusted_offset, dna_sequence_ptr->length(),
                         del_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;
  }

  ContigOffset_t sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);

  // Check the reference.
  if (del_ptr->reference() != dna_sequence_ptr->at(sequence_offset)) {

    ExecEnv::log().warn("mutateSingleDelete(), reference base: {} does not match sequence base: {} at sequence offset: {}",
                        DNA5::convertToChar(del_ptr->reference()),
                        DNA5::convertToChar(dna_sequence_ptr->at(sequence_offset)),
                        sequence_offset);

  }

// Important - This routine is paired with adjustIndelOffsets(), which should be called BEFORE calling (this routine) updateIndelAccounting()
// Important - update the indel offset accounting structure AFTER the actual indel offset has been calculated with adjustIndelOffsets().// adjust for any prior indels.
  if (not updateIndelAccounting(del_ptr, indel_accounting_map)) {

    return false;

  }

  // Mutate the sequence
  dna_sequence_ptr->deleteSubSequence(sequence_offset, del_ptr->size());

  return true;

}


// Mutate the DNA sequence using Insert variants
bool kgl::VariantMutation::mutateInserts(const OffsetVariantMap& snp_variant_map,
                                         ContigOffset_t contig_offset,
                                         std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                                         IndelAccountingMap& indel_accounting_map) {

  return true;

}
