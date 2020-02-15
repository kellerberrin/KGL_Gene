//
// Created by kellerberrin on 30/08/18.
//

#ifndef KGL_VARIANT_FACTORY_VCF_CIGAR_H
#define KGL_VARIANT_FACTORY_VCF_CIGAR_H



#include "kgl_utility.h"
#include "kgl_variant_evidence.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_vcf_parse_impl.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF parser. Low-level implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class ParseCigar : public ParseVCFImpl {

public:

  ParseCigar(std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                 std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                 const std::string& vcf_file_name) : ParseVCFImpl(vcf_population_ptr, genome_db_ptr, vcf_file_name) {}
  ~ParseCigar() override = default;

  bool parseCigarItems(const std::string& genome_name,
                       std::shared_ptr<const ContigFeatures> contig_ptr,
                       const std::vector<CigarEditItem>& parsed_cigar,
                       ContigOffset_t contig_offset,
                       const std::string& reference,
                       const std::string& alternate,
                       std::shared_ptr<const VariantEvidence> evidence_ptr,
                       size_t& record_variants);

private:


};


}   // end namespace



#endif //KGL_KGL_VARIANT_FACTORY_VCF_CIGAR_H
