//
// Created by kellerberrin on 22/01/18.
//

#ifndef KGL_VARIANT_FACTORY_GATKVCF_IMPL_H
#define KGL_VARIANT_FACTORY_GATKVCF_IMPL_H



#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_single.h"
#include "kgl_variant_factory_readvcf_impl.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <seqan/vcf_io.h>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class GATKVCFImpl : public ParseVCFImpl {

public:

  GATKVCFImpl(const std::string& genome_name,
              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
              const std::string& vcf_file_name,
              Phred_t variant_quality) :   genome_name_(genome_name),
                                           genome_db_ptr_(genome_db_ptr),
                                           vcf_file_name_(vcf_file_name),
                                           variant_quality_(variant_quality) {

    reader_ptr_ = std::make_shared<VCFReaderMT<GATKVCFImpl>>(vcf_file_name, this, &GATKVCFImpl::ProcessVCFRecord);
    genome_single_variants_ = GenomeVariant::emptyGenomeVariant(genome_name_, genome_db_ptr_);

  }

  ~GATKVCFImpl() override = default;

  std::shared_ptr<GenomeVariant> readParseGATKVcfFile();

  void ProcessVCFRecord(const seqan::VcfRecord& record_ptr);

private:

  std::shared_ptr<GenomeVariant> processParseGATKVcfFile(const std::string& genome_name,
                                                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                         const std::string& vcf_file_name,
                                                         Phred_t variant_quality);

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

  bool parseDelete(const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& reference,
                   const std::string& alternate,
                   ContigOffset_t contig_offset,
                   size_t& variant_count) const;

  const std::string& genome_name_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  const std::string& vcf_file_name_;
  Phred_t variant_quality_;
  std::shared_ptr<GenomeVariant> genome_single_variants_;

  std::shared_ptr<VCFReaderMT<GATKVCFImpl>> reader_ptr_;


};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_GATKVCF_IMPL_H
