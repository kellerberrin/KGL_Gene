//
// Created by kellerberrin on 22/01/18.
//

#ifndef KGL_VARIANT_FACTORY_GATKVCF_IMPL_H
#define KGL_VARIANT_FACTORY_GATKVCF_IMPL_H



#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_single.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <seqan/vcf_io.h>


namespace kgl = kellerberrin::genome;
namespace bt = boost;

class kgl::VcfFactory::GATKVCFImpl : public kgl::VcfFactory::ParseVCFImpl {

public:

  GATKVCFImpl() = default;
  ~GATKVCFImpl() override = default;

  std::shared_ptr<GenomeVariant> readParseGATKVcfFile(const std::string& genome_name,
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
                      size_t& record_variants) const;

  bool parseSNP(const std::string& variant_source,
                std::shared_ptr<const ContigFeatures> contig_ptr,
                std::shared_ptr<GenomeVariant> genome_variants,
                Phred_t quality,
                const std::string& info,
                const std::string& reference,
                const std::string& alternate,
                ContigOffset_t contig_offset,
                size_t& variant_count) const;

  bool parseInsert(const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& reference,
                   const std::string& alternate,
                   ContigOffset_t contig_offset,
                   size_t& variant_count) const;


};






#endif //KGL_VARIANT_FACTORY_GATKVCF_IMPL_H