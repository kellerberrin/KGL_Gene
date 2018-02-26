//
// Created by kellerberrin on 25/02/18.
//

#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_pf3k_impl.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <seqan/vcf_io.h>


namespace kgl = kellerberrin::genome;
namespace bt = boost;



bool kgl::Pf3kVCFImpl::readParsePf3kVariants() {

  ExecEnv::log().info("readParsePf3kVariants(), parsing file: {}", vcf_file_name_);

  VCFReaderMT<Pf3kVCFImpl> reader(vcf_file_name_, this, &Pf3kVCFImpl::ProcessVCFRecord);

  reader.readVCFFile();

  return true;

}

// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::ProcessVCFRecord(const seqan::VcfRecord& record_ptr)
{

  ++vcf_variant_count_;

  if (vcf_variant_count_ % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("VCF file, generated: {} variants", vcf_variant_count_);

  }

}
