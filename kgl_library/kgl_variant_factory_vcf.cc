//
// Created by kellerberrin on 26/12/17.
//

#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_fbvcf_impl.h"
#include "kgl_variant_factory_gatkvcf_impl.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto VcfFactory::FreeBayesVCFImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::VcfFactory::VcfFactory() : fb_vcf_impl_ptr_(std::make_unique<kgl::VcfFactory::FreeBayesVCFImpl>()),
                                gatk_vcf_impl_ptr_(std::make_unique<kgl::VcfFactory::GATKVCFImpl>()) {}
kgl::VcfFactory::~VcfFactory() {}  // DO NOT DELETE or USE DEFAULT. Required here because of incomplete pimpl types.


std::shared_ptr<kgl::GenomeVariant> kgl::VcfFactory::readParseFreeBayesVcf(const std::string &genome_name,
                                                                           std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                                           const std::string &vcf_file_name,
                                                                           Phred_t variant_quality) const {

  return fb_vcf_impl_ptr_->readParseFreeBayesVcfFile(genome_name,
                                                     genome_db_ptr,
                                                     vcf_file_name,
                                                     variant_quality);

}


std::shared_ptr<kgl::GenomeVariant> kgl::VcfFactory::readParseGATKVcf(const std::string &genome_name,
                                                                      std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                                      const std::string &vcf_file_name,
                                                                      Phred_t variant_quality) const {

  return gatk_vcf_impl_ptr_->readParseGATKVcfFile(genome_name,
                                                  genome_db_ptr,
                                                  vcf_file_name,
                                                  variant_quality);

}

