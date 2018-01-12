//
// Created by kellerberrin on 26/12/17.
//

#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_vcf_impl.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto VcfFactory::VcfFileImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::VcfFactory::VcfFactory() : vcf_file_impl_ptr_(std::make_unique<kgl::VcfFactory::VcfFileImpl>()) {}
kgl::VcfFactory::~VcfFactory() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.


std::shared_ptr<kgl::GenomeVariant> kgl::VcfFactory::readParseVcf(const std::string& genome_name,
                                                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                                  const std::string& vcf_file_name,
                                                                  Phred_t variant_quality) const {

  return vcf_file_impl_ptr_->readParseVcfFile(genome_name,
                                              genome_db_ptr,
                                              vcf_file_name,
                                              variant_quality);

}

