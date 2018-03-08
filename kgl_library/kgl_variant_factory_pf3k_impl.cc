//
// Created by kellerberrin on 25/02/18.
//

#include "kgl_ploidy_analysis.h"
#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_pf3k_impl.h"


namespace kgl = kellerberrin::genome;

// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::ProcessVCFRecord(const seqan::VcfRecord& vcf_record) {

  try {

    TryVCFRecord(vcf_record);

  }
  catch(const std::exception& e) {

    ExecEnv::log().error("ProcessVCFRecord(), Exception: {} thrown record ignored", e.what());

  }


}

// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::TryVCFRecord(const seqan::VcfRecord& vcf_record) {

  ++vcf_variant_count_;

  std::shared_ptr<PloidyAnalysis> ploidy_ptr = std::dynamic_pointer_cast<PloidyAnalysis>(pop_variant_ptr_);

  if (ploidy_ptr) {

    process_ploidy_.PloidyVCFRecord(vcf_record, ploidy_ptr, contigId(vcf_record.rID));

  }

  if (vcf_variant_count_ % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("VCF file, processed: {} variants", vcf_variant_count_);

  }

}
