//
// Created by kellerberrin on 16/10/17.
//

#include <sstream>
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;

std::string kgl::ReadCountFilter::filterName() const {

  std::ostringstream oss;
  oss << "Minimum Read Count >= " << read_count_;
  return oss.str();

}


std::string kgl::MutantProportionFilter::filterName() const {

  std::ostringstream oss;
  oss << "Minimum Read Proportion >= " << mutant_proportion_;
  return oss.str();

}


bool kgl::InCDSFilter::implementFilter(const Variant& variant) const {

  std::shared_ptr<ContigFeatures> contig_ptr;
  if (genome_db_ptr_->getContigSequence(variant.contigId(), contig_ptr)) {

    std::vector<std::shared_ptr<CDSFeature>> cds_ptr_vec;
    return contig_ptr->findOffsetCDS(variant.contigOffset(), cds_ptr_vec);

  } else {

    ExecEnv::log().error("Variant contig: {} not found in genome database", variant.contigId());

  }

  return false;

}


std::string kgl::InCDSFilter::filterName() const {

  return "Variant in CDS";

}
