//
// Created by kellerberrin on 11/03/18.
//

#ifndef KGL_VARIANT_FACTORY_VCF_CIGAR_IMPL_H
#define KGL_VARIANT_FACTORY_VCF_CIGAR_IMPL_H




#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_vcf_parse_impl.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF (freebayes) parser. Low-level implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class ParseCigarImpl : public ParseVCFImpl {

public:

  ParseCigarImpl(std::shared_ptr<VCFPopulation> vcf_population_ptr,
                 std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                 const std::string& vcf_file_name,
                 Phred_t variant_quality) : ParseVCFImpl(vcf_population_ptr, genome_db_ptr, vcf_file_name, variant_quality) {}
  ~ParseCigarImpl() override = default;

  bool parseCigarItems(const std::string& genome_name,
                       std::shared_ptr<const ContigFeatures> contig_ptr,
                       const std::vector<CigarEditItem>& parsed_cigar,
                       ContigOffset_t contig_offset,
                       const std::string& reference,
                       const std::string& alternate,
                       Phred_t quality,
                       const std::string& info,
                       size_t& record_variants);

private:


  // Parse 1M ... XM in the cigar.
  bool parseCheck(size_t cigar_count,
                  std::shared_ptr<const ContigFeatures> contig_ptr,
                  const std::string& reference,
                  const std::string& alternate,
                  size_t& reference_index,
                  size_t& alternate_index,
                  ContigOffset_t& contig_offset) const;

  // Parse 1X ... XX in the cigar.
  bool parseSNP(size_t cigar_count,
                const std::string& variant_source,
                std::shared_ptr<const ContigFeatures> contig_ptr,
                Phred_t quality,
                const std::string& info,
                const std::string& reference,
                const std::string& alternate,
                size_t& reference_index,
                size_t& alternate_index,
                ContigOffset_t& contig_offset,
                size_t& variant_count);

  // Parse 1I ... XI in the cigar.
  bool parseInsert(size_t cigar_count,
                   const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& alternate,
                   ContigOffset_t contig_offset,
                   size_t& alternate_index,
                   size_t& variant_count);

  // Parse 1D ... XD in the cigar.
  bool parseDelete(size_t cigar_count,
                   const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& reference,
                   size_t& reference_index,
                   ContigOffset_t& contig_offset,
                   size_t& variant_count);




};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_FACTORY_VCF_CIGAR_IMPL_H
