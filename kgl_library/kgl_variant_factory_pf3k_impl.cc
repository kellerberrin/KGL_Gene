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



// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::ProcessVCFRecord(const seqan::VcfRecord& record_ptr)
{

  ++vcf_variant_count_;

  if (vcf_variant_count_ % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("VCF file, generated: {} variants", vcf_variant_count_);

  }

  if (getGenomeNames().size() != seqan::length(record_ptr.genotypeInfos)) {

    ExecEnv::log().warn("Genome Name Size: {}, Genotype count: {}",
                        getGenomeNames().size(), seqan::length(record_ptr.genotypeInfos));

  }


  static bool first_record = true;

  if (first_record) {

    first_record = false;

    std::cout << "Alt:" << record_ptr.alt << std::endl;
    std::cout << "Pos:" << record_ptr.beginPos << std::endl;
    std::cout << "Filter:" << record_ptr.filter << std::endl;
    std::cout << "Format:" << record_ptr.format << std::endl;
    std::cout << "Id:" << record_ptr.id << std::endl;
    std::cout << "Info:" << record_ptr.info << std::endl;
    std::cout << "Qual:" << record_ptr.qual << std::endl;
    std::cout << "Ref:" << record_ptr.ref << std::endl;
    std::cout << "Rid:" << record_ptr.rID << std::endl;

    auto it = seqan::begin(record_ptr.genotypeInfos);
    auto end = seqan::end(record_ptr.genotypeInfos);

    while (it != end)
    {

      std::cout << "Genotype:" << *it << std::endl;
      ++it;

    }


  }


}
