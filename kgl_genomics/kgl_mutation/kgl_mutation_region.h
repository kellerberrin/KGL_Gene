//
// Created by kellerberrin on 2/08/23.
//

#ifndef KGL_MODIFIED_DNA_REGION_H
#define KGL_MODIFIED_DNA_REGION_H



#include "kgl_variant_db.h"
#include "kgl_mutation_db.h"
#include "kgl_mutation_variant_offset.h"

#include <map>
#include <memory>
#include <vector>
#include <sstream>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A variant map modifies a DNA region.
//
// This object is supplied with a reference contig and a map of canonical variants over a specific contig region.
// The object creates a modified sequence specified by the variant map and also produces an offset map that
// details how indels have modified the original offset structure of the modified sequence.
// This permits any feature such as an exon or intron to be mapped onto the modified sequence.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ModifiedDNARegion {

public:

  ModifiedDNARegion(const RegionVariantMap &variant_map, const std::shared_ptr<const ContigReference>& contig_ptr) {

    mutateDNA(variant_map, contig_ptr);

  }
  ~ModifiedDNARegion() = default;

  [[nodiscard]] const VariantMutationOffset& mutationOffsetMap() const { return variant_modification_offset_; }
  [[nodiscard]] const DNA5SequenceLinear& mutatedSequence() const { return modified_sequence_; }

private:

  VariantMutationOffset variant_modification_offset_;
  DNA5SequenceLinear modified_sequence_;

  // These routines perform the actual mutation of the DNA sequence.
  // Modification variants, SNP and indel, in the map must be in canonical format.
  void mutateDNA(const RegionVariantMap &variant_map, const std::shared_ptr<const ContigReference>& contig_ptr);

  // Modify the sequence with SNP and indel variants. An indel will return a non-zero signed offset. A +ve
  // signed offset is returned if an insert variant, a -ve offset is returned if a delete variant.
  // Note that delete indels can precede the offset of the modified region and can overlap the end of the of
  // the modified region.
  [[nodiscard]] std::pair<SignedOffset_t, bool> modifyVariant(const std::shared_ptr<const Variant>& variant_ptr,
                                                              SignedOffset_t offset_adjust);

  // Mutate a sequence by adding and subtracting subsequences at the designated adjusted offset
  [[nodiscard]] std::pair<SignedOffset_t, bool> performModification(const std::shared_ptr<const Variant>& variant_ptr,
                                                                    ContigOffset_t offset_adjust);

  // Called for a delete variant that precedes the region but the deleted region overlaps the modified region.
  [[nodiscard]] std::pair<SignedOffset_t, bool> precedingVariant(const std::shared_ptr<const Variant>& variant_ptr,
                                                                 SignedOffset_t adjusted_offset);

};



}   // end namespace


#endif //KGL_MODIFIED_DNA_REGION_H
