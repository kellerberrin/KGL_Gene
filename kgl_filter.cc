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
// Filter SNPs to coding sequences only.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::InCDSFilter::implementFilter(const Variant& variant) const {

  return variant.type() == VariantGenomeType::CDS_CODING;

}


std::string kgl::InCDSFilter::filterName() const {

  return "Variant in CDS";

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::ContigFilter::implementFilter(const Variant& variant) const {

    return variant.contigId() == contig_ident_;

}

std::string kgl::ContigFilter::filterName() const {

  return "Variant in Contig: " + contig_ident_;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular gene.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GeneFilter::implementFilter(const Variant& variant) const {

  std::shared_ptr<ContigFeatures> contig_ptr;
  if (genome_db_ptr_->getContigSequence(variant.contigId(), contig_ptr)) {

    std::vector<std::shared_ptr<CDSFeature>> cds_ptr_vec;
    if (not contig_ptr->findOffsetCDS(variant.contigOffset(), cds_ptr_vec)) {

      return false;

    } else {

      std::shared_ptr<kgl::Feature> gene_ptr = cds_ptr_vec.front()->getGene();
      if (not gene_ptr) {

        ExecEnv::log().error("Variant GeneFilter; Gene not found for CDS: {}", cds_ptr_vec.front()->id());
        return false;

      } else {

        return gene_ptr->id() == gene_ident_;

      }

    }


  } else {

    ExecEnv::log().error("Variant contig: {} not found in genome database", variant.contigId());

  }

  return false;

}


std::string kgl::GeneFilter::filterName() const {

  return "Variant in Gene: " + gene_ident_;

}

