//
// Created by kellerberrin on 25/02/18.
//

#include "kgl_variant_factory_record_vcf_impl.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_pf3k_impl.h"


namespace kgl = kellerberrin::genome;



// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::ProcessVCFRecord(const seqan::VcfRecord& vcf_record) {

  static bool first_record = true;

  ++vcf_variant_count_;

  ParseVCFRecord recordParser(vcf_record, contigId(vcf_record.rID), genome_db_ptr_);

  if (vcf_variant_count_ % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("VCF file, generated: {} variants", vcf_variant_count_);
    first_record = true;

  }

  if (getGenomeNames().size() != seqan::length(vcf_record.genotypeInfos)) {

    ExecEnv::log().warn("Genome Name Size: {}, Genotype count: {}",
                        getGenomeNames().size(), seqan::length(vcf_record.genotypeInfos));

  }


  if (first_record) {


    std::cout << "Alt:" << vcf_record.alt << std::endl;
    std::cout << "Pos:" << vcf_record.beginPos << std::endl;
    std::cout << "Filter:" << vcf_record.filter << std::endl;
    std::cout << "Format:" << vcf_record.format << std::endl;
    std::cout << "Id:" << vcf_record.id << std::endl;
    std::cout << "Info:" << vcf_record.info << std::endl;
    std::cout << "Qual:" << vcf_record.qual << std::endl;
    std::cout << "Ref:" << vcf_record.ref << std::endl;
    std::cout << "Rid:" << vcf_record.rID << std::endl;

/*
    auto it = seqan::begin(vcf_record.genotypeInfos);
    auto end = seqan::end(vcf_record.genotypeInfos);

    while (it != end)
    {

      std::cout << "Genotype:" << *it << std::endl;
      ++it;

    }

*/

    auto it = seqan::begin(vcf_record.genotypeInfos);
    auto end = seqan::end(vcf_record.genotypeInfos);
    size_t genotype_count = 0;
    size_t valid_genotype_count = 0;

    while (it != end)
    {

      ParseVCFGenotype genotype_parser(*it);

      if (genotype_parser.formatCount() == recordParser.requiredFormatSize()) {

        if (genotype_parser.getFormatChar(recordParser.PLOffset(), *it) != PL_CHECK_ZERO_) {

          ++valid_genotype_count;

          ExecEnv::log().info("Alleles : {}, GT : {}, GQ : {}, PL : {}",
                              recordParser.alleles().size(),
                              genotype_parser.getFormatString(recordParser.GTOffset(), *it),
                              genotype_parser.getFormatString(recordParser.GQOffset(), *it),
                              genotype_parser.getFormatString(recordParser.PLOffset(), *it));

        }

      }

      ++genotype_count;
      ++it;

    }

    ExecEnv::log().info("Genotype count: {}, Valid Genotype count: {}", genotype_count, valid_genotype_count);

  }

  first_record = false;

}
