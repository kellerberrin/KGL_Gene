//
// Created by kellerberrin on 30/08/18.
//

#include "kgl_variant_vcf.h"
#include "kgl_variant_factory_vcf_cigar.h"


namespace kgl = kellerberrin::genome;


bool kgl::ParseCigar::parseCigarItems(const std::string& genome_name,
                                      std::shared_ptr<const ContigFeatures> contig_ptr,
                                      const std::vector<CigarEditItem>&, // parsed_cigar,  // Cigar not used
                                      ContigOffset_t contig_offset,
                                      const std::string& reference_text,
                                      const std::string& alternate_text,
                                      std::shared_ptr<const VariantEvidence> evidence_ptr,
                                      size_t& record_variants)  {

  StringDNA5 reference_str(reference_text);
  StringDNA5 alternate_str(alternate_text);

  std::shared_ptr<const Variant> variant_ptr(std::make_shared<VCFVariant>(genome_name,
                                                                          contig_ptr->contigId(),
                                                                          VariantSequence::UNPHASED,
                                                                          contig_offset,
                                                                          evidence_ptr,
                                                                          std::move(reference_str),
                                                                          std::move(alternate_str)));

  record_variants = addThreadSafeGenomeVariant(variant_ptr);

  return true;

}

