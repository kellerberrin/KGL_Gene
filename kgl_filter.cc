//
// Created by kellerberrin on 16/10/17.
//

#include <sstream>
#include <memory>
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set the minimum read count SNP generation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string kgl::ReadCountFilter::filterName() const {

  std::ostringstream oss;
  oss << "Minimum Read Count >= " << read_count_;
  return oss.str();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set the minimum mutant read proportion in a candidate SNP.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string kgl::MutantProportionFilter::filterName() const {

  std::ostringstream oss;
  oss << "Minimum Read Proportion >= " << mutant_proportion_;
  return oss.str();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variantss to coding sequences only.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::InCDSFilter::implementFilter(const Variant& variant) const {

  return variant.type() == VariantSequenceType::CDS_CODING;

}


std::string kgl::InCDSFilter::filterName() const {

  return "Variant in CDS";

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::ContigFilter::implementFilter(const Variant& variant) const {

    return variant.contigId() == contig_ident_;

}

std::string kgl::ContigFilter::filterName() const {

  return "Variant in Contig: " + contig_ident_;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a particular gene.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GeneFilter::implementFilter(const Variant& variant) const {

  GeneVector gene_ptr_vec;
  if (variant.type() != VariantSequenceType::NON_CODING) {

    for (const auto& gene_ptr : variant.geneMembership()) {

      if (gene_ptr->id() == gene_ident_) {

        return true;

      }

    }

  }

  return false;

}


std::string kgl::GeneFilter::filterName() const {

  return "Variant in Gene: " + gene_ident_;

}

