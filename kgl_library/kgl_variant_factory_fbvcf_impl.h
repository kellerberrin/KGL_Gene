//
// Created by kellerberrin on 22/01/18.
//

#ifndef KGL_VARIANT_FACTORY_FBVCF_IMPL_H
#define KGL_VARIANT_FACTORY_FBVCF_IMPL_H




#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_single.h"

#include <seqan/vcf_io.h>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF (freebayes) parser. Low-level implementation. Do not include this file in any source files except the following:
// kgl_variant_factory_vcf.cc
// kgl_variant_factory_vcf_fbimpl.cc
// kgl_variant_factory_vcf_impl.cc
// VcfFactory::FreeBayesVCFImpl does all the heavy lifting using the 3rd party libraries, seqan and boost.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class FreeBayesVCFImpl : public ParseVCFImpl {

public:

  FreeBayesVCFImpl() = default;
  ~FreeBayesVCFImpl() override = default;

  std::shared_ptr<GenomeVariant> readParseFreeBayesVcfFile(const std::string& genome_name,
                                                           std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                           const std::string& vcf_file_name,
                                                           Phred_t variant_quality);

private:


  bool parseVcfRecord(const std::string& genome_name,
                      const seqan::VcfRecord& record,
                      std::shared_ptr<const ContigFeatures> contig_ptr,
                      std::shared_ptr<GenomeVariant> genome_variants,
                      Phred_t variant_quality,
                      bool& quality_ok,
                      size_t& variant_count) const;

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
                std::shared_ptr<GenomeVariant> genome_variants,
                Phred_t quality,
                const std::string& info,
                const std::string& reference,
                const std::string& alternate,
                size_t& reference_index,
                size_t& alternate_index,
                ContigOffset_t& contig_offset,
                size_t& variant_count) const;

  // Parse 1I ... XI in the cigar.
  bool parseInsert(size_t cigar_count,
                   const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& alternate,
                   ContigOffset_t contig_offset,
                   size_t& alternate_index,
                   size_t& variant_count) const;

  // Parse 1D ... XD in the cigar.
  bool parseDelete(size_t cigar_count,
                   const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& reference,
                   size_t& reference_index,
                   ContigOffset_t& contig_offset,
                   size_t& variant_count) const;


};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_FBVCF_IMPL_H
