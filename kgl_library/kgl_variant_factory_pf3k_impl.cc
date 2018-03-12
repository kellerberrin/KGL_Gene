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

  } else {

    ParseRecord(vcf_record, contigId(vcf_record.rID));

  }

  if (vcf_variant_count_ % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("VCF file, processed: {} variants", vcf_variant_count_);

  }

}


// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::ParseRecord(const seqan::VcfRecord& vcf_record, const ContigId_t& contig_id) {

  ParseVCFRecord recordParser(vcf_record, contig_id, genome_db_ptr_); //Each vcf record.

  VCFInfoField info_key_value_map(seqan::toCString(vcf_record.info));  // Each vcf record.

  if (getGenomeNames().size() != seqan::length(vcf_record.genotypeInfos)) {

    ExecEnv::log().warn("Genome Name Size: {}, Genotype count: {}",
                        getGenomeNames().size(), seqan::length(vcf_record.genotypeInfos));

  }

  auto begin = seqan::begin(vcf_record.genotypeInfos);
  auto end = seqan::end(vcf_record.genotypeInfos);
  size_t genotype_count = 0;

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

        std::string DP_text = genotype_parser.getFormatString(recordParser.DPOffset(), *it);
        if (DP_text.find_first_not_of(DIGITS_) !=  std::string::npos) {

          continue;

        }

        size_t DP_value = std::stoll(DP_text);
        size_t A_allele = std::stoll(gt_vector[0]);
        size_t B_allele = std::stoll(gt_vector[1]);

        const std::string &genome_name = getGenomeNames()[genotype_count];
        bool valid_record = vqslod >= MIN_VQSLOD_QUALITY_
                            and recordParser.quality() >= MIN_QUALITY_
                            and GQ_value >= MIN_GQ_QUALITY_
                            and DP_value >= MIN_DEPTH_;


        if (valid_record and A_allele > 0) {

          const std::string allele = recordParser.alleles()[A_allele-1];
          if (allele != UPSTREAM_ALLELE_) {

            // process A allele
            std::vector<CigarEditItem> parsed_cigar;
            ParseVCFMiscImpl::generateEditVector(recordParser.reference(), allele, parsed_cigar);
            std::string info;
            size_t record_variants;
            std::shared_ptr<GenomeVariant> genome_variant;
            if (not thread_safe_population_.getCreateGenomeVariant(genome_name, genome_db_ptr_, genome_variant)) {

              ExecEnv::log().error("Could not find or create genome: {}", genome_name);
              continue;

            }

            parseCigarItems(genome_name,
                            recordParser.contigPtr(),
                            genome_variant,
                            parsed_cigar,
                            recordParser.offset(),
                            recordParser.reference(),
                            allele,
                            recordParser.quality(),
                            info,
                            record_variants);

            variant_count_ += parsed_cigar.size();

          }

        }

        if (valid_record and B_allele > 0) {

          const std::string allele = recordParser.alleles()[B_allele-1];
          if (allele != UPSTREAM_ALLELE_) {

            // process B allele
            std::vector<CigarEditItem> parsed_cigar;
            ParseVCFMiscImpl::generateEditVector(recordParser.reference(), allele, parsed_cigar);
            std::string info;
            size_t record_variants;
            std::shared_ptr<GenomeVariant> genome_variant;
            if (not thread_safe_population_.getCreateGenomeVariant(genome_name, genome_db_ptr_, genome_variant)) {

              ExecEnv::log().error("Could not find or create genome: {}", genome_name);
              continue;

            }

            parseCigarItems(genome_name,
                            recordParser.contigPtr(),
                            genome_variant,
                            parsed_cigar,
                            recordParser.offset(),
                            recordParser.reference(),
                            allele,
                            recordParser.quality(),
                            info,
                            record_variants);

            variant_count_ += parsed_cigar.size();

          }

        }

      }

    }
    ++record_count_;
    ++genotype_count;

    if (record_count_ % 1000000 == 0) {

      ExecEnv::log().info("Processed :{} records, variants: {}, variants: {}",
                          static_cast<size_t>(record_count_), variant_count_, thread_safe_population_.variantCount());

      for (auto allele : recordParser.alleles()) {

        ExecEnv::log().info("Reference: {}, Alternate: {}, Cigar: {}",
                            recordParser.reference(), allele, ParseVCFMiscImpl::generateCigar(recordParser.reference(), allele));

      }

    }

  }

}
