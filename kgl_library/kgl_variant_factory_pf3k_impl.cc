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

  static bool first_record = true;

  ++vcf_variant_count_;

  std::shared_ptr<PloidyAnalysis> ploidy_ptr = std::dynamic_pointer_cast<PloidyAnalysis>(pop_variant_ptr_);

  ParseVCFRecord recordParser(vcf_record, contigId(vcf_record.rID), genome_db_ptr_);


  if (vcf_variant_count_ % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("VCF file, generated: {} variants, Quality Failed: {}, Vqslod Failed: {}, GQ Failed: {}",
                        vcf_variant_count_, static_cast<size_t>(quality_failed_), static_cast<size_t>(vqslod_failed_), static_cast<size_t>(GQ_failed_));

  }

  if (recordParser.quality() < variant_quality_) {

    ++quality_failed_;
    return;

  }

  if (getGenomeNames().size() != seqan::length(vcf_record.genotypeInfos)) {

    ExecEnv::log().warn("Genome Name Size: {}, Genotype count: {}",
                        getGenomeNames().size(), seqan::length(vcf_record.genotypeInfos));

  }

  VCFInfoField info_key_value_map(seqan::toCString(vcf_record.info));


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

          ExecEnv::log().info("Diploid Genotypes : {}", diploid_genotypes_.genotypeText(recordParser.alleles().size()));
          ExecEnv::log().info("Alleles : {}, GT : {}, AD: {}, DP: {}, GQ : {}, PL : {}",
                              recordParser.alleles().size(),
                              genotype_parser.getFormatString(recordParser.GTOffset(), *it),
                              genotype_parser.getFormatString(recordParser.ADOffset(), *it),
                              genotype_parser.getFormatString(recordParser.DPOffset(), *it),
                              genotype_parser.getFormatString(recordParser.GQOffset(), *it),
                              genotype_parser.getFormatString(recordParser.PLOffset(), *it));

        }

      }

      ++genotype_count;
      ++it;

    }

    ExecEnv::log().info("Genotype count: {}, Valid Genotype count: {}", genotype_count, valid_genotype_count);

    first_record = false;

  } // First record

  if (ploidy_ptr) {

    auto it = seqan::begin(vcf_record.genotypeInfos);
    auto end = seqan::end(vcf_record.genotypeInfos);
    size_t genotype_count = 0;
    bool one_error = false;

    while (it != end)
    {

      ParseVCFGenotype genotype_parser(*it);

      if (genotype_parser.formatCount() == recordParser.requiredFormatSize()) {

        if (genotype_parser.getFormatChar(recordParser.PLOffset(), *it) != PL_CHECK_ZERO_) {

          std::string vqslod_text;
          if ( not info_key_value_map.getInfoField(VQSLOD_INFO_FIELD_, vqslod_text)) {

            ExecEnv::log().error("Parsing Pf3k VCF, Info field : {} not found", VQSLOD_INFO_FIELD_);
            continue;

          }

          double vqslod = std::atof(vqslod_text.c_str());

          if (vqslod < MIN_VQSLOD_QUALITY_) {

            ++vqslod_failed_;
            return;

          }


          std::vector<std::string> gt_vector;
          ParseVCFMiscImpl::tokenize(genotype_parser.getFormatString(recordParser.GTOffset(), *it), GT_FIELD_SEPARATOR_, gt_vector);

          if (gt_vector.size() != 2) {

            ExecEnv::log().error("Parsing Pf3k VCF, PT format field: {} is not diploid.",
                                 genotype_parser.getFormatString(recordParser.GTOffset(), *it));
            continue;

          }


          double GQ_value = std::stof(genotype_parser.getFormatString(recordParser.GQOffset(), *it));
          if (GQ_value < MIN_GQ_QUALITY_) {

            ++GQ_failed_;
            return;

          }

          size_t A_allele = std::stoll(gt_vector[0]);
          size_t B_allele = std::stoll(gt_vector[1]);

          std::vector<std::string> pl_vector;
          ParseVCFMiscImpl::tokenize(genotype_parser.getFormatString(recordParser.PLOffset(), *it), PL_FIELD_SEPARATOR_, pl_vector);

          size_t genotype_permutations = diploid_genotypes_.generateGenotype(recordParser.alleles().size()).size();
          if (pl_vector.size() != genotype_permutations) {

            ExecEnv::log().error("Parsing Pf3k VCF, PL format field: {} permutation mismatch",
                                 genotype_parser.getFormatString(recordParser.PLOffset(), *it));
            continue;

          }

          std::vector<size_t> pl_phred_vec;
          for (auto pl_item : pl_vector) {

            pl_phred_vec.push_back(std::stoll(pl_item));

          }

          bool found_zero = false;
          size_t pl_index = 0;
          size_t A_pl_allele;
          size_t B_pl_allele;
          for (auto pl_pred : pl_phred_vec) {

            if (pl_pred == 0) {

              A_pl_allele = diploid_genotypes_.generateGenotype(recordParser.alleles().size())[pl_index].first;
              B_pl_allele = diploid_genotypes_.generateGenotype(recordParser.alleles().size())[pl_index].second;
              found_zero = true;

            }

            ++pl_index;

          }

          if (found_zero) {


            if ((A_allele != A_pl_allele or B_allele != B_pl_allele) and not one_error and GQ_value > 0) {

              AutoMutex error_info_mutex(mutex_);

              ExecEnv::log().error("Parsing Pf3k VCF Mismatching alleles, A allele: {}, A PL Allele : {}, B allele : {}, B PL Allele : {}",
                                   A_allele, A_pl_allele, B_allele, B_pl_allele);
              ExecEnv::log().info("Diploid Genotypes : {}", diploid_genotypes_.genotypeText(recordParser.alleles().size()));
              ExecEnv::log().info("Ref (0): {}, Alleles(1,...): {}",  seqan::toCString(vcf_record.ref), seqan::toCString(vcf_record.alt));
              ExecEnv::log().info("Quality: {}, Vqslod: {}, Alleles: {}; GT :{}; AD:{}; DP:{}; GQ :{}; PL :{}",
                                  vcf_record.qual,
                                  vqslod,
                                  recordParser.alleles().size(),
                                  genotype_parser.getFormatString(recordParser.GTOffset(), *it),
                                  genotype_parser.getFormatString(recordParser.ADOffset(), *it),
                                  genotype_parser.getFormatString(recordParser.DPOffset(), *it),
                                  genotype_parser.getFormatString(recordParser.GQOffset(), *it),
                                  genotype_parser.getFormatString(recordParser.PLOffset(), *it));
              one_error = true;
              continue;

            }

          } else {

            ExecEnv::log().error("Parsing Pf3k VCF, Did not find zero PL field: {}",
                                 genotype_parser.getFormatString(recordParser.PLOffset(), *it));

            continue;

          }

          const std::string &genome = getGenomeNames()[genotype_count];
          bool homozygous = (A_allele == B_allele);
          bool hq_homozygous = homozygous and GQ_value >= HQ_GQ_PLOIDY_;
          bool heterozygous = (A_allele != B_allele);
          bool hq_heterozygous = heterozygous and GQ_value >= HQ_GQ_PLOIDY_;
          if (not ploidy_ptr->addPloidyRecord(genome, homozygous, hq_homozygous, heterozygous, hq_heterozygous)) {

            ExecEnv::log().error("Cannot add ploidy record");

          }

        }

      }
      ++ploidy_count_;
      ++genotype_count;
      ++it;

      if (ploidy_count_ % 1000000 == 0) {

        ExecEnv::log().info("Processed :{} ploidy records", ploidy_count_);

      }

    }

  }


}
