//
// Created by kellerberrin on 25/02/18.
//

#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_pf3k_impl.h"


namespace kgl = kellerberrin::genome;


// Process VCF header information.
void kgl::Pf3kVCFImpl::processVCFHeader(const VcfHeaderInfo& header_info) {

  // Investigate header.
  ActiveContigMap active_contig_map;
  if (not ParseVCFMiscImpl::parseVcfHeader(genome_db_ptr_, header_info, active_contig_map, false)) {

    ExecEnv::log().error("Problem parsing header information in VCF file: {}. No variants processed.", vcf_file_name_);

  }

}

void kgl::Pf3kVCFImpl::readParseVCFImpl() {



  // multi-threaded
  readVCFFile();
  // single threaded



}

// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) {

  try {

    TryVCFRecord(vcf_record_count, vcf_record);

  }
  catch(const std::exception& e) {

    ExecEnv::log().error("ProcessVCFRecord(), Exception: {} thrown record ignored", e.what());

  }


}

// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::TryVCFRecord(size_t vcf_record_count, const VcfRecord& record) {

  ++vcf_variant_count_;

  ParseRecord(vcf_record_count, record);

  if (vcf_variant_count_ % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("VCF file, processed: {} variants", vcf_variant_count_);

  }

}


// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::ParseRecord(size_t vcf_record_count, const VcfRecord& record) {

  ParseVCFRecord recordParser(record, record.contig_id, genome_db_ptr_); //Each vcf record.

  VCFInfoField info_key_value_map(record.info);  // Each vcf record.

  if (getGenomeNames().size() != record.genotypeInfos.size()) {

    ExecEnv::log().warn("Genome Name Size: {}, Genotype count: {}",
                        getGenomeNames().size(), record.genotypeInfos.size());

  }

  size_t genotype_count = 0;

  for (auto const& genotype : record.genotypeInfos)
  {

    ParseVCFGenotype genotype_parser(genotype);

    if (genotype_parser.formatCount() == recordParser.requiredFormatSize()) {

      if (genotype_parser.getFormatChar(recordParser.PLOffset(), genotype) != PL_CHECK_ZERO_
          and genotype_parser.getFormatChar(recordParser.PLOffset(), genotype) != PL_CHECK_DOT_) {

        std::string vqslod_text;
        bool vqslod_found = true;
        if ( not info_key_value_map.getInfoField(VQSLOD_INFO_FIELD_, vqslod_text)) {

          vqslod_found = false;

        }

        double vqslod = 0.0;
        if (vqslod_found) {

          if (vqslod_text.find_first_not_of(FLOAT_DIGITS_) !=  std::string::npos) {

            ExecEnv::log().info("Non-numeric vqslod text: {}", vqslod_text);
            continue;

          }
          vqslod = std::atof(vqslod_text.c_str());

        }

        std::vector<std::string> gt_vector;

        if (not ParseVCFMiscImpl::tokenize(genotype_parser.getFormatString(recordParser.GTOffset(), genotype), GT_FIELD_SEPARATOR_, gt_vector)) {

          ExecEnv::log().error("Parsing Pf3k VCF, Problem tokenizing PT format field: {}.",
                               genotype_parser.getFormatString(recordParser.GTOffset(), genotype));
          continue;

        }

        if (gt_vector.size() != 2) {

          ExecEnv::log().error("Parsing Pf3k VCF, PT format field: {} is not diploid.",
                               genotype_parser.getFormatString(recordParser.GTOffset(), genotype));
          continue;

        }

        std::string GQ_text = genotype_parser.getFormatString(recordParser.GQOffset(), genotype);
        if (GQ_text.find_first_not_of(DIGITS_) !=  std::string::npos) {

          ExecEnv::log().error("Non-numeric GQ_text: {}", GQ_text);
          continue;

        }
        double GQ_value = std::stof(GQ_text);

        std::string DP_text = genotype_parser.getFormatString(recordParser.DPOffset(), genotype);
        if (DP_text.find_first_not_of(DIGITS_) !=  std::string::npos) {

          continue;

        }

        size_t DP_value = std::stoll(DP_text);
        size_t A_allele = std::stoll(gt_vector[0]);
        size_t B_allele = std::stoll(gt_vector[1]);

        // Get ad allele depths.
        std::vector<std::string> ad_vector;

        if (not ParseVCFMiscImpl::tokenize(genotype_parser.getFormatString(recordParser.ADOffset(), genotype), AD_FIELD_SEPARATOR_, ad_vector)) {

          ExecEnv::log().error("Parsing Pf3k VCF, Problem tokenizing AD format field");
          continue;

        }

        // Allele depths should be the number of alleles + the reference
        if (ad_vector.size() != (recordParser.alleles().size() + 1)) {

          ExecEnv::log().error("Parsing Pf3k VCF, Expected:{} AD allele depths, actually:{}",
                               (recordParser.alleles().size() + 1), ad_vector.size());
          continue;

        }

        std::vector<size_t> ad_count_vector;
        size_t ad_total_count = 0;
        for (auto const& depth_count_text : ad_vector) {

          if (depth_count_text.find_first_not_of(DIGITS_) !=  std::string::npos) {

            ExecEnv::log().error("Parsing Pf3k VCF, Expected numeric value AD allele depths, actually:{}", depth_count_text);
            continue;

          }

          size_t ad_count = std::stoll(depth_count_text);
          ad_total_count += ad_count;
          ad_count_vector.push_back(ad_count);

        }

        const std::string &genome_name = getGenomeNames()[genotype_count];

        bool valid_record = recordParser.quality() >= MIN_QUALITY_
                            and GQ_value >= MIN_GQ_QUALITY_
                            and DP_value >= MIN_DEPTH_;

        if (vqslod_found) {

          valid_record = valid_record and vqslod >= MIN_VQSLOD_QUALITY_;

        }

        if (valid_record and A_allele > 0) {

          size_t allele_index = A_allele - 1;   // 0 is the reference

          const std::string allele = recordParser.alleles()[allele_index];

          size_t ref_count = ad_count_vector[0];
          size_t alt_count = ad_count_vector[A_allele];

          // VCF variants with zero alt+ref counts are flagged as spanning (downstream)
          // deletion. There is a spanning upstream delete. The downstream variant is ignored.
          bool downstream_variant = ref_count == 0 and alt_count == 0;

          if (allele != UPSTREAM_ALLELE_ and not downstream_variant) {

            // Evidence object
            std::shared_ptr<VariantEvidence> evidence_ptr(std::make_shared<CountEvidence>(ref_count,
                                                                                          alt_count,
                                                                                          DP_value,
                                                                                          GQ_value,
                                                                                          recordParser.quality(),
                                                                                          vcf_record_count));
            // process A allele
            std::vector<CigarEditItem> parsed_cigar = ParseVCFMiscImpl::generateEditVector(recordParser.reference(), allele);
            size_t record_variants;

            if (not parseCigarItems( genome_name,
                                     recordParser.contigPtr(),
                                     parsed_cigar,
                                     recordParser.offset(),
                                     recordParser.reference(),
                                     allele,
                                     evidence_ptr,
                                  record_variants)) {

              ExecEnv::log().error("Parsing Pf3k VCF, Problem parsing A allele CIGAR items");

            }

            variant_count_ += parsed_cigar.size();

          }

        }

        if (valid_record and B_allele > 0) {

          size_t allele_index = B_allele - 1;   // 0 is the reference

          const std::string allele = recordParser.alleles()[allele_index];

          size_t ref_count = ad_count_vector[0];
          size_t alt_count = ad_count_vector[B_allele];

          // VCF variants with zero alt+ref counts are flagged as spanning (downstream)
          // deletion. There is a spanning upstream delete. The downstream variant is ignored.
          bool downstream_variant = ref_count == 0 and alt_count == 0;

          if (allele != UPSTREAM_ALLELE_ and not downstream_variant) {

            // Evidence object
            std::shared_ptr<VariantEvidence> evidence_ptr(std::make_shared<CountEvidence>(ref_count,
                                                                                          alt_count,
                                                                                          DP_value,
                                                                                          GQ_value,
                                                                                          recordParser.quality(),
                                                                                          vcf_record_count));

            // process B allele
            std::vector<CigarEditItem> parsed_cigar = ParseVCFMiscImpl::generateEditVector(recordParser.reference(), allele);
            size_t record_variants;

            if (not parseCigarItems( genome_name,
                                     recordParser.contigPtr(),
                                     parsed_cigar,
                                     recordParser.offset(),
                                     recordParser.reference(),
                                     allele,
                                     evidence_ptr,
                                  record_variants)) {

              ExecEnv::log().error("Parsing Pf3k VCF, Problem parsing B allele CIGAR items");

            }

            variant_count_ += parsed_cigar.size();

          }

        }

      }

    }
    ++record_count_;
    ++genotype_count;

    if (record_count_ % 1000000 == 0) {

      ExecEnv::log().info("Processed :{} records, total variants: {}", static_cast<size_t>(record_count_), variant_count_);
      ExecEnv::log().info("Contig: {}, offset: {}", recordParser.contigPtr()->contigId(), recordParser.offset());

    }

  }

}
