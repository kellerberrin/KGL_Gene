//
// Created by kellerberrin on 13/10/17.
//

#include <ostream>
#include "kgl_variant.h"
#include "kgl_patterns.h"
#include "kgl_filter.h"
#include "kgl_variant_db.h"
#include "kgl_sequence_offset.h"


namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


const kgl::PhaseId_t kgl::VariantSequence::UNPHASED;


std::string kgl::VariantSequence::genomeOutput(char delimiter, VariantOutputIndex output_index) const {

  std:: stringstream ss;
// Contig.
  ss << genomeId() << delimiter
     << contigId() << delimiter;
  if (phaseId() == UNPHASED) {

    ss << "UNPHASED" << delimiter;

  } else {

    ss << "PHASE:" << static_cast<size_t>(phaseId()) << delimiter;

  }
  ss << offsetOutput(offset(), output_index) << delimiter;

  return ss.str();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string kgl::Variant::name() const {

 switch(variantType()) {

   case VariantType::SNP: return "SNP";

   case VariantType::DELETE: return "DEL";

   case VariantType::INSERT: return "INS";

   case VariantType::COMPOUND_DELETE: return "CDEL";

   case VariantType::COMPOUND_INSERT: return "CINS";

 }

 return "NOT_IMPLEMENTED";  // Not reached, to keep the compiler happy.

}


// Relevant for compound variants.
// Used to see if variants are mutating the same section of DNA.
bool kgl::Variant::offsetOverlap(const Variant& cmp_var) const {

  // First check if the variants are on the same contig.
  if (contigId() != cmp_var.contigId()) {

    return false;

  }

  // On the same contig so check for overlap.
  ContigOffset_t start = offset();
  ContigOffset_t end = offset() + size() - 1;

  ContigOffset_t cmp_start = cmp_var.offset();
  ContigOffset_t cmp_end = cmp_var.offset() + size() - 1;

  bool overlap_offset = ((start <= cmp_start) and (cmp_start <= end)) or ((start <= cmp_end) and (cmp_end <= end));

  // The overlap response depends on type of variant overlapping.
  // If Insert / Insert - abnormal complain and split
  // If Insert / Delete - normal no split.
  // If Insert / SNP - normal do not split.
  // If Delete / Delete - abnormal complain and split
  // If Delete and SNP - abnormal complain and split.
  // If SNO and SNP - normal split.
  if (overlap_offset) {

    if (isInsert()) {

      if (cmp_var.isInsert()) {

        // complain if both have the same offset
        if (offset() == cmp_var.offset()) {

          ExecEnv::log().warn("offsetOverlap(), Two insert variants with the same offset: {} overlaps Insert variant: {}",
                              output(' ', VariantOutputIndex::START_0_BASED, true),
                              cmp_var.output(' ', VariantOutputIndex::START_0_BASED, true));
          return true; // Complain and Split
        }
        return false; // OK and no split

      } else if (cmp_var.isDelete()) {

        return false; // OK and no split

      } else if (cmp_var.isSNP()) {

        return false; // OK and no split

      } else {

        ExecEnv::log().warn("offsetOverlap(), Unexpected variant type: {}",
                             cmp_var.output(' ', VariantOutputIndex::START_0_BASED, true));
        return false; // no split

      }

    } else if (isDelete()) {

      if (cmp_var.isInsert()) {

        // complain
        ExecEnv::log().warn("offsetOverlap(), Delete variant: {} overlaps Insert variant: {}",
                             output(' ', VariantOutputIndex::START_0_BASED, true),
                             cmp_var.output(' ', VariantOutputIndex::START_0_BASED, true));
        return true; // and split

      } else if (cmp_var.isDelete()) {

        // complain
        ExecEnv::log().warn("offsetOverlap(), Delete variant: {} overlaps Delete variant: {}",
                             output(' ', VariantOutputIndex::START_0_BASED, true),
                             cmp_var.output(' ', VariantOutputIndex::START_0_BASED, true));
        return true; // and split

      } else if (cmp_var.isSNP()) {

        // complain
        ExecEnv::log().warn("offsetOverlap(), Delete variant: {} overlaps SNP variant: {}",
                             output(' ', VariantOutputIndex::START_0_BASED, true),
                             cmp_var.output(' ', VariantOutputIndex::START_0_BASED, true));
        return true; // and split

      } else {

        ExecEnv::log().error("offsetOverlap(), Unexpected variant type: {}",
                             cmp_var.output(' ', VariantOutputIndex::START_0_BASED, true));
        return false; // no split

      }


    } else if (isSNP()) {

      if (cmp_var.isInsert()) {

        return false; // OK and no split

      } else if (cmp_var.isDelete()) {

        // complain
        ExecEnv::log().warn("offsetOverlap(), SNP variant: {} overlaps Delete variant: {}",
                            output(' ', VariantOutputIndex::START_0_BASED, true),
                            cmp_var.output(' ', VariantOutputIndex::START_0_BASED, true));
        return true; // and split

      } else if (cmp_var.isSNP()) {

        return true; // OK and split

      } else {

        ExecEnv::log().error("offsetOverlap(), Unexpected variant type: {}",
                             cmp_var.output(' ', VariantOutputIndex::START_0_BASED, true));
        return false; // no split

      }

    } else {

      ExecEnv::log().error("offsetOverlap(), Unexpected variant type: {}",
                           output(' ', VariantOutputIndex::START_0_BASED, true));
      return false; // no split

    }

  }

  return false; // no overlap.

}