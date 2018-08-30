//
// Created by kellerberrin on 30/08/18.
//


#include "kgl_variant_factory_vcf_cigar.h"
#include "kgl_variant_factory_compound.h"


namespace kgl = kellerberrin::genome;


bool kgl::ParseCigar::parseCigarItems(const std::string& genome_name,
                                      std::shared_ptr<const ContigFeatures> contig_ptr,
                                      const std::vector<CigarEditItem>& parsed_cigar,
                                      ContigOffset_t begin_contig_offset,
                                      const std::string& reference,
                                      const std::string& alternate,
                                      std::shared_ptr<const VariantEvidence> evidence_ptr,
                                      size_t& record_variants)  {




  return true;

}

