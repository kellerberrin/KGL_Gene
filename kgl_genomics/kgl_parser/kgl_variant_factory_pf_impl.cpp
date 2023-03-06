//
// Created by kellerberrin on 25/02/18.
//

#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_pf_impl.h"
#include "kgl_variant.h"


namespace kgl = kellerberrin::genome;


// Process VCF header information.
void kgl::PfVCFImpl::processVCFHeader(const VCFHeaderInfo& header_info) {

  // Investigate header.
  VCFContigMap vcf_contig_map;
  VCFInfoRecordMap vcf_info_map;
  if (not VCFParseHeader::parseVcfHeader( header_info, vcf_contig_map, vcf_info_map)) {

    ExecEnv::log().error("PfVCFImpl::processVCFHeader, Problem parsing header information in VCF file. No variants processed.");

  }

  // Store all the available VCF info fields.
  evidence_factory_.availableInfoFields(vcf_info_map);

  if (VCFParseHeader::checkVCFReferenceContigs(vcf_contig_map, genome_db_ptr_)) {

    ExecEnv::log().info("VCF File and Reference Genome Contig Sizes All Match");

  } else {

    ExecEnv::log().info("VCF File and Reference Genome Contig Size/Name Mis-Match");

  }

  // Pre-initializes the UnphasedPopulation with a list of genomes/contigs.
  // Since some genomes may not have a variant.
  // The initialization of the reader above generates a list of genomes in the VCF file.
  setupPopulationStructure(genome_db_ptr_);

}

void kgl::PfVCFImpl::readParseVCFImpl(const std::string &vcf_file_name) {

  // multi-threaded
  readVCFFile(vcf_file_name);
  // single threaded

}

// This is multithreaded code called from the reader defined above.
void kgl::PfVCFImpl::ProcessVCFRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr) {

  try {

    ParseRecord(std::move(vcf_record_ptr));

  }
  catch(const std::exception& e) {

    ExecEnv::log().error("PfVCFImpl::ProcessVCFRecord, Unexpected Exception: {} thrown record ignored", e.what());

  }

}


// This is multithreaded code called from the reader defined above.
void kgl::PfVCFImpl::ParseRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr) {

  ParseVCFRecord recordParser(*vcf_record_ptr, genome_db_ptr_); //Each vcf record.
  if (not recordParser.parseResult()) {

    ExecEnv::log().warn("PfVCFImpl::ParseRecord; Problem parsing VCF record");
    return;

  }

  // Parse the info fields into a map.
  // For performance reasons the info field is std::moved - don't reference again.
  auto mutable_info = const_cast<std::string&>(vcf_record_ptr->info);
  std::shared_ptr<const DataMemoryBlock> info_evidence_ptr = evidence_factory_.createVariantEvidence(std::move(mutable_info));  // Each vcf record.

  if (getGenomeNames().size() != vcf_record_ptr->genotypeInfos.size()) {

    ExecEnv::log().warn("PfVCFImpl::ParseRecord; Genome Name Size: {}, Genotype count: {}",
                        getGenomeNames().size(), vcf_record_ptr->genotypeInfos.size());

  }

  auto GT_offset_opt = recordParser.formatIndex(ParseVCFRecord::FORMAT_GT);
  auto GQ_offset_opt = recordParser.formatIndex(ParseVCFRecord::FORMAT_GQ);
  auto DP_offset_opt = recordParser.formatIndex(ParseVCFRecord::FORMAT_DP);
  auto AD_offset_opt = recordParser.formatIndex(ParseVCFRecord::FORMAT_AD);
  // Require GT and AD format fields.
  if (not GT_offset_opt or not AD_offset_opt) {

    ExecEnv::log().error("PfVCFImpl::ParseRecord; format: {} does not contain required format fields: {} and {}",
                         vcf_record_ptr->format, ParseVCFRecord::FORMAT_GT, ParseVCFRecord::FORMAT_AD);
    return;
  }

  size_t genotype_count = 0;
  for (auto const& genotype : vcf_record_ptr->genotypeInfos)
  {

    std::vector<std::string_view> genotype_formats = Utility::viewTokenizer(genotype, FORMAT_SEPARATOR_);

    // Require GT format field.
    if (genotype_formats.size() <= GT_offset_opt.value()) {

      ExecEnv::log().error( "PfVCFImpl::ParseRecord; record: {}, genotype format: {} does not contain required format field: {}",
                            vcf_record_ptr->line_number, genotype, ParseVCFRecord::FORMAT_GT);
      continue;

    }

    std::string GT_format(genotype_formats[GT_offset_opt.value()]);
    std::vector<std::string> gt_vector = Utility::charTokenizer(GT_format, GT_FIELD_SEPARATOR_CHAR_);
    if (gt_vector.size() != DIPLOID_) {

      ExecEnv::log().error("PfVCFImpl::ParseRecord; GT format field: {} is not diploid.", GT_format, genotype);
      continue;

    }

    size_t A_allele{0};
    if (gt_vector[0].find_first_not_of(DIGITS_) ==  std::string::npos) {

      A_allele = std::stoll(gt_vector[0]);

    }

    size_t B_allele{0};
    if (gt_vector[0].find_first_not_of(DIGITS_) ==  std::string::npos) {

      B_allele = std::stoll(gt_vector[1]);

    }

    // If there are no alt alleles then skip.
    if (A_allele == 0 and B_allele == 0) {

      continue;

    }

    // Get the GQ field if it exists.
    double GQ_value{0.0};
    if (GQ_offset_opt) {

      if (GQ_offset_opt.value() < genotype_formats.size()) {

        std::string GQ_text(genotype_formats[GQ_offset_opt.value()]);

        if (GQ_text != MISSING_VALUE_) {

          if (GQ_text.find_first_not_of(FLOAT_DIGITS_) !=  std::string::npos) {

            ExecEnv::log().error("PfVCFImpl::ParseRecord; Non-numeric GQ_text: {}", GQ_text);

          } else {

            GQ_value = std::stof(GQ_text);

          }

        }

      }

    }

    // Get the DP field if it exists.
    size_t DP_value{0};
    if (DP_offset_opt) {

      if (DP_offset_opt.value() < genotype_formats.size()) {

        std::string DP_text(genotype_formats[DP_offset_opt.value()]);

        if (DP_text != MISSING_VALUE_) {

          if (DP_text.find_first_not_of(DIGITS_) !=  std::string::npos) {

            ExecEnv::log().error("PfVCFImpl::ParseRecord; Non-numeric DP_text: {}", DP_text);

          } else {

            DP_value = std::stoll(DP_text);

          }

        }

      }

    }

    // Require AD format field.
    if (genotype_formats.size()  <= AD_offset_opt.value()) {

      ExecEnv::log().error( "PfVCFImpl::ParseRecord; record: {}, genotype format: {} does not contain required format field: {}",
                            vcf_record_ptr->line_number, genotype, ParseVCFRecord::FORMAT_AD);
      continue;

    }

    std::vector<size_t> ad_count_vector;
    if (AD_offset_opt) {

      // Get ad allele depths.
      std::string AD_text(genotype_formats[AD_offset_opt.value()]);
      std::vector<std::string> ad_vector = Utility::charTokenizer(AD_text, AD_FIELD_SEPARATOR_CHAR_);

      // Allele depths should be the number of alleles + the reference
      if (ad_vector.size() != (recordParser.alleles().size() + 1)) {

        ExecEnv::log().error("PfVCFImpl::ParseRecord; Expected:{} AD allele depths, actually:{}",
                             (recordParser.alleles().size() + 1), ad_vector.size());
        continue;

      }

      size_t ad_total_count = 0;
      for (auto const& depth_count_text : ad_vector) {

        if (depth_count_text.find_first_not_of(DIGITS_) !=  std::string::npos) {

          ExecEnv::log().error("PfVCFImpl::ParseRecord; Expected numeric value AD allele depths, actually:{}", depth_count_text);
          continue;

        }

        size_t ad_count = std::stoll(depth_count_text);
        ad_total_count += ad_count;
        ad_count_vector.push_back(ad_count);

      }

    } else {

      ExecEnv::log().error("PfVCFImpl::ParseRecord; AD format field not found");
      continue;

    }

    const std::string &genome_name = getGenomeNames()[genotype_count];

    if (A_allele > 0) {

      uint32_t allele_count = recordParser.alleles().size();
      uint32_t allele_index = A_allele - 1;   // 0 is the reference

      const std::string allele = recordParser.alleles()[allele_index];

      size_t ref_count = ad_count_vector[0];
      size_t alt_count = ad_count_vector[A_allele];

      // VCF variants with zero alt+ref counts are flagged as spanning (downstream)
      // deletion. There is a spanning upstream delete. The downstream variant is ignored.
      bool downstream_variant = ref_count == 0 and alt_count == 0;

      if (allele != UPSTREAM_ALLELE_ and not downstream_variant) {

        // Format Evidence object
        std::shared_ptr<FormatData> format_data_ptr(std::make_shared<FormatData>(ref_count,
                                                                                 alt_count,
                                                                                 DP_value,
                                                                                 GQ_value,
                                                                                 recordParser.quality()));
        // Setup the evidence object
        VariantEvidence evidence(vcf_record_ptr->line_number,
                                 unphased_population_ptr_->dataSource(),
                                 recordParser.passedFilter(),
                                 info_evidence_ptr,
                                 format_data_ptr,
                                 allele_index,
                                 allele_count);

        if (not createAddVariant(genome_name,
                                 recordParser.contigPtr(),
                                 recordParser.offset(),
                                 vcf_record_ptr->id,
                                 recordParser.reference(),
                                 allele,
                                 evidence)) {

          ExecEnv::log().error("PfVCFImpl::ParseRecord; Problem parsing A allele");

        }

        ++variant_count_;

      }

    }

    if (B_allele > 0) {

      uint32_t allele_index = B_allele - 1;   // 0 is the reference
      uint32_t allele_count = recordParser.alleles().size();
      const std::string allele = recordParser.alleles()[allele_index];

      size_t ref_count = ad_count_vector[0];
      size_t alt_count = ad_count_vector[B_allele];

      // VCF variants with zero alt+ref counts are flagged as spanning (downstream)
      // deletion. There is a spanning upstream delete. The downstream variant is ignored.
      bool downstream_variant = ref_count == 0 and alt_count == 0;

      if (allele != UPSTREAM_ALLELE_ and not downstream_variant) {

        // Format Evidence object
        std::shared_ptr<FormatData> format_data_ptr(std::make_shared<FormatData>(ref_count,
                                                                                 alt_count,
                                                                                 DP_value,
                                                                                 GQ_value,
                                                                                 recordParser.quality()));
        // Setup the evidence object.
        VariantEvidence evidence(vcf_record_ptr->line_number,
                                 unphased_population_ptr_->dataSource(),
                                 recordParser.passedFilter(),
                                 info_evidence_ptr,
                                 format_data_ptr,
                                 allele_index,
                                 allele_count);

        if (not createAddVariant(genome_name,
                                 recordParser.contigPtr(),
                                 recordParser.offset(),
                                 vcf_record_ptr->id,
                                 recordParser.reference(),
                                 allele,
                                 evidence)) {

          ExecEnv::log().error("PfVCFImpl::ParseRecord; Problem parsing B allele");

        }

        ++variant_count_;

      }

    }

    // Next genome name.
    ++genotype_count;

  }

  if (vcf_record_ptr->line_number % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("Processed :{} records, total variants: {}", vcf_record_ptr->line_number, variant_count_);
    ExecEnv::log().info("Contig: {}, offset: {}", recordParser.contigPtr()->contigId(), recordParser.offset());

  }

}


// Set up the genomes/contigs first rather than on-the-fly.
// Some genomes may have no variants (e.g. the model/reference genome 3D7)
// and thus these genomes/contigs would not be created on-the-fly.
void kgl::PfVCFImpl::setupPopulationStructure(const std::shared_ptr<const GenomeReference>& genome_db_ptr) {

  ExecEnv::log().info("PfVCFImpl::setupPopulationStructure; Creating a population of {} genomes and {} contigs",
                      getGenomeNames().size(), genome_db_ptr->getMap().size());

  for (auto const& genome_id : getGenomeNames())  {

    std::optional<std::shared_ptr<GenomeDB>> genome_opt = unphased_population_ptr_->getCreateGenome(genome_id);
    if (not genome_opt) {
      // Terminate runtime.
      ExecEnv::log().critical("PfVCFImpl::setupPopulationStructure; Could not create genome: {} in the unphased population", genome_id);

    }

    for (auto const& contig : genome_db_ptr->getMap()) {

      if (not genome_opt.value()->getCreateContig(contig.first)) {
        // Terminate runtime.
        ExecEnv::log().critical("PfVCFImpl::setupPopulationStructure; Could not create contig: {} in genome: {} in the unphased population", contig.first, genome_id);

      }

    }

  }

}

bool kgl::PfVCFImpl::createAddVariant(const std::string& genome_name,
                                      const std::shared_ptr<const ContigReference>& contig_ptr,
                                      ContigOffset_t contig_offset,
                                      const std::string& identifier,
                                      const std::string& reference_text,
                                      const std::string& alternate_text,
                                      const VariantEvidence& evidence)  {

  StringDNA5 reference_str(reference_text);
  StringDNA5 alternate_str(alternate_text);

  std::shared_ptr<const Variant> variant_ptr(std::make_shared<const Variant>( contig_ptr->contigId(),
                                                                              contig_offset,
                                                                              VariantPhase::UNPHASED,
                                                                              identifier,
                                                                              std::move(reference_str),
                                                                              std::move(alternate_str),
                                                                              evidence));

  return addThreadSafeVariant(variant_ptr, genome_name);

}


bool kgl::PfVCFImpl::addThreadSafeVariant(const std::shared_ptr<const Variant>& variant_ptr, const GenomeId_t& genome) const {

  std::vector<GenomeId_t> genome_vector;
  genome_vector.push_back(genome);

  // The thread mutex is within the addVariant() function.
  return unphased_population_ptr_->addVariant(variant_ptr, genome_vector);

}
