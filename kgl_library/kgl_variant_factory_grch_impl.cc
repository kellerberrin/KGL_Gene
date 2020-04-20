//
// Created by kellerberrin on 20/4/20.
//

#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_grch_impl.h"


namespace kgl = kellerberrin::genome;


void kgl::GrchVCFImpl::processVCFHeader(const VcfHeaderInfo& header_info) {


}


void kgl::GrchVCFImpl::readParseVCFImpl() {



  // multi-threaded
  readVCFFile();
  // single threaded

  for (auto const& [contig_id, count] : contig_count_) {

    ExecEnv::log().info("GrchVCFImpl; Contig id: {}, length: {}, variant count :{}", contig_id, count.first, count.second);

  }

}


void kgl::GrchVCFImpl::ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) {

  contig_count_[vcf_record.contig_id].second++;

}

