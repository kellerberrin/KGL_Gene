//
// Created by kellerberrin on 7/03/18.
//


#include "kgl_ploidy_vcf_Pf3k_impl.h"
#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_pf3k_impl.h"
#include "kgl_sequence_compare_impl.h"


namespace kgl = kellerberrin::genome;



// This is multithreaded code called from the reader defined above.
void kgl::ProcessPloidy::PloidyVCFRecord(const seqan::VcfRecord& vcf_record,
                                         std::shared_ptr<PloidyAnalysis> ploidy_ptr,
                                         const ContigId_t& contig_id) {

  ParseVCFRecord recordParser(vcf_record, contig_id, genome_db_ptr_); //Each vcf record.

  VCFInfoField info_key_value_map(seqan::toCString(vcf_record.info));  // Each vcf record.

  if (genome_names_.size() != seqan::length(vcf_record.genotypeInfos)) {

    ExecEnv::log().warn("Genome Name Size: {}, Genotype count: {}",
                        genome_names_.size(), seqan::length(vcf_record.genotypeInfos));

  }

  for (auto allele : recordParser.alleles()) {

    std::vector<SequenceEditType> edit_vector;
    SequenceComparison().generateEditVector(recordParser.reference(), allele, edit_vector);

  }

  auto begin = seqan::begin(vcf_record.genotypeInfos);
  auto end = seqan::end(vcf_record.genotypeInfos);
  size_t genotype_count = 0;
  bool one_error = false;

  for (auto it = begin; it !=end; ++it)
  {

    ParseVCFGenotype genotype_parser(*it);

    if (genotype_parser.formatCount() == recordParser.requiredFormatSize()) {

      if (genotype_parser.getFormatChar(recordParser.PLOffset(), *it) != PL_CHECK_ZERO_
          and genotype_parser.getFormatChar(recordParser.PLOffset(), *it) != PL_CHECK_DOT_) {

        std::string vqslod_text;
        if ( not info_key_value_map.getInfoField(VQSLOD_INFO_FIELD_, vqslod_text)) {

          ExecEnv::log().error("Parsing Pf3k VCF, Info field : {} not found", VQSLOD_INFO_FIELD_);
          continue;

        }

        if (vqslod_text.find_first_not_of(FLOAT_DIGITS_) !=  std::string::npos) {

          ExecEnv::log().info("Non-numeric vqslod text: {}", vqslod_text);
          continue;

        }
        double vqslod = std::atof(vqslod_text.c_str());

//           if (vqslod < MIN_VQSLOD_QUALITY_) {
//
//            ++vqslod_failed_;
//            return;
//
//           }

        std::vector<std::string> gt_vector;
        ParseVCFMiscImpl::tokenize(genotype_parser.getFormatString(recordParser.GTOffset(), *it), GT_FIELD_SEPARATOR_, gt_vector);

        if (gt_vector.size() != 2) {

          ExecEnv::log().error("Parsing Pf3k VCF, PT format field: {} is not diploid.",
                               genotype_parser.getFormatString(recordParser.GTOffset(), *it));
          continue;

        }

        std::string GQ_text = genotype_parser.getFormatString(recordParser.GQOffset(), *it);
        if (GQ_text.find_first_not_of(DIGITS_) !=  std::string::npos) {

          ExecEnv::log().info("Non-numeric GQ_text: {}", GQ_text);
          continue;

        }
        double GQ_value = std::stof(GQ_text);
//          if (GQ_value < MIN_GQ_QUALITY_) {
//
//            ++GQ_failed_;
//            return;
//
//          }
//
        std::string DP_text = genotype_parser.getFormatString(recordParser.DPOffset(), *it);
        if (DP_text.find_first_not_of(DIGITS_) !=  std::string::npos) {

//          ExecEnv::log().info("Non-numeric DP_text: {}", DP_text);
          continue;

        }
        size_t DP_value = std::stoll(DP_text);
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


        const std::string &genome = genome_names_[genotype_count];
        bool homozygous = (A_allele == B_allele);
        bool hq_homozygous = homozygous
                             and GQ_value >= HQ_GQ_PLOIDY_
                             and recordParser.quality() >= MIN_QUALITY_
                             and GQ_value >= MIN_GQ_QUALITY_
                             and DP_value >= MIN_DEPTH_
                             and recordParser.isSNP();
        bool heterozygous = (A_allele != B_allele);
        bool hq_heterozygous = heterozygous
                               and GQ_value >= HQ_GQ_PLOIDY_
                               and recordParser.quality() >= MIN_QUALITY_
                               and GQ_value >= MIN_GQ_QUALITY_
                               and DP_value >= MIN_DEPTH_
                               and recordParser.isSNP();


        double ratio_count = 0;
        if (hq_homozygous or hq_heterozygous) {

          std::vector<std::string> ad_vector;
          ParseVCFMiscImpl::tokenize(genotype_parser.getFormatString(recordParser.ADOffset(), *it), AD_FIELD_SEPARATOR_, ad_vector);

          if (ad_vector.size() != recordParser.alleles().size() + 1) {

            ExecEnv::log().error("Parsing Pf3k VCF, AD format field: {} size: {} not equal size alleles + 1",
                                 genotype_parser.getFormatString(recordParser.ADOffset(), *it), recordParser.alleles().size() + 1);
            continue;

          }

          if (A_allele >= ad_vector.size() or B_allele >= ad_vector.size()) {

            ExecEnv::log().error("A allele index: {}, B allele index: {} out of range for AD vector: {}",
                                 A_allele, B_allele, genotype_parser.getFormatString(recordParser.ADOffset(), *it));
            continue;

          }

          if (A_allele != B_allele and A_allele > 0 and B_allele > 0) {

            ++different_alleles_;

          }


          double A_proportion = std::stof(ad_vector[A_allele]);
          double B_proportion = std::stof(ad_vector[B_allele]);

          double max_count = A_proportion >= B_proportion ? A_proportion : B_proportion;
          double min_count = A_proportion <= B_proportion ? A_proportion : B_proportion;

          if (max_count <= 0.0) {

            ExecEnv::log().info("Max count for AD vector: {} is zero, DP value: {}; ref: {}, alt: {}",
                                genotype_parser.getFormatString(recordParser.ADOffset(), *it), DP_value,
                                seqan::toCString(vcf_record.ref), seqan::toCString(vcf_record.alt));

            ExecEnv::log().info("Genotype record: {}", seqan::toCString(*it));

            ratio_count = 0.0;

          } else {

            ratio_count = min_count / max_count;

          }

        }

        if (not ploidy_ptr->addPloidyRecord(genome, homozygous, hq_homozygous, heterozygous, hq_heterozygous, ratio_count)) {

          ExecEnv::log().error("Cannot add ploidy record");

        }

      }

    }
    ++ploidy_count_;
    ++genotype_count;

    if (ploidy_count_ % 1000000 == 0) {

      ExecEnv::log().info("Processed :{} ploidy records, different alleles: {}", ploidy_count_, different_alleles_);

      for (auto allele : recordParser.alleles()) {

        std::vector<SequenceEditType> edit_vector;
        SequenceComparison().generateEditVector(recordParser.reference(), allele, edit_vector);
        std::string edit_string;
        for (auto edit : edit_vector) {

          switch(edit) {

            case SequenceEditType::UNCHANGED:
              edit_string += "M,";
              break;

            case SequenceEditType::INSERT:
              edit_string += "I,";
              break;

            case SequenceEditType::DELETE:
              edit_string += "D,";
              break;

            case SequenceEditType::CHANGED:
              edit_string += "X,";
              break;

          }


        }
        ExecEnv::log().info("Reference: {}, Alternate: {}, Cigar: {}, Edit string: {}",
                            recordParser.reference(), allele, SequenceComparison().generateCigar(recordParser.reference(), allele), edit_string);

      }

    }

  }

}
