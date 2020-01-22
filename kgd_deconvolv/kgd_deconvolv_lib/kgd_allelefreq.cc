//
// Created by kellerberrin on 19/05/18.
//
#include "kgd_allelefreq.h"

namespace kgd = kellerberrin::deconvolv;


// Get or Create ContigAlleles
kgd::ContigAlleles& kgd::GenomeAlleles::getCreateContigAllele(const ContigId_t& contig_id) {

  auto result = contig_allele_map_.find(contig_id);

  if (result != contig_allele_map_.end()) {

    return result->second;

  }

  auto insert_result = contig_allele_map_.insert(ContigAlleleMap::value_type(contig_id, contig_id));

  if (not insert_result.second) {

    // Something very wrong with memory alloc.
    DeconvolvApp::log().critical("GenomeAlleles.getCreateContigAllele(); could not add ContigAlleles:{}", contig_id);

  }

  return insert_result.first->second;

}


void kgd::GenomeAlleles::addAlleleFreq(const ContigId_t& contig_id, ContigOffset_t offset, AlleleFreq_t frequency) {


  ContigAlleles& contig_allele = getCreateContigAllele(contig_id);

  contig_allele.addAlleleFreq(AlleleOffset(offset, frequency));

}
