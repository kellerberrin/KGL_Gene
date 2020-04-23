//
// Created by kellerberrin on 20/4/20.
//

#include "kgl_variant_factory_vcf_parse_impl.h"
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

  // Parse the info fields into an array and assign to a shared ptr.
//  std::shared_ptr<VCFInfoField> info_key_value_map_ptr(std::make_shared<VCFInfoField>(vcf_record.info));  // Each vcf record.
  std::scoped_lock<std::mutex> lock(mutex_);

  auto mutable_info = const_cast<std::string&>(vcf_record.info);
  std::shared_ptr<std::string> info_ptr = std::make_shared<std::string>(std::move(mutable_info));

  contig_count_[vcf_record.contig_id].second++;

}

