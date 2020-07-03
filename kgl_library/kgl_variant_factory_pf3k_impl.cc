//
// Created by kellerberrin on 25/02/18.
//

#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_pf3k_impl.h"
#include "kgl_variant.h"


namespace kgl = kellerberrin::genome;


// Process VCF header information.
void kgl::Pf3kVCFImpl::processVCFHeader(const VcfHeaderInfo& header_info) {

  // Investigate header.
  VCFContigMap vcf_contig_map;
  VCFInfoRecordMap vcf_info_map;
  if (not VCFParseHeader::parseVcfHeader( header_info, vcf_contig_map, vcf_info_map)) {

    ExecEnv::log().error("Pf3kVCFImpl::processVCFHeader, Problem parsing header information in VCF file. No variants processed.");

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

void kgl::Pf3kVCFImpl::readParseVCFImpl() {

  // multi-threaded
  readVCFFile();
  // single threaded

}

// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) {

  try {

    ParseRecord(vcf_record_count, vcf_record);

  }
  catch(const std::exception& e) {

    ExecEnv::log().error("Pf3kVCFImpl::ProcessVCFRecord, Unexpected Exception: {} thrown record ignored", e.what());

  }

}


// This is multithreaded code called from the reader defined above.
void kgl::Pf3kVCFImpl::ParseRecord(size_t vcf_record_count, const VcfRecord& record) {

  ParseVCFRecord recordParser(record); //Each vcf record.

  if (not recordParser.parseRecord(record.contig_id, genome_db_ptr_)) {

    ExecEnv::log().warn("Pf3kVCFImpl::ParseRecord, Problem parsing VCF record");
    return;

  }

  // Parse the info fields into a map.
  // For performance reasons the info field is std::moved - don't reference again.
  auto mutable_info = const_cast<std::string&>(record.info);
  InfoDataEvidence info_evidence_opt = evidence_factory_.createVariantEvidence(std::move(mutable_info));  // Each vcf record.

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

        std::vector<std::string> gt_vector = Utility::char_tokenizer(genotype_parser.getFormatString(recordParser.GTOffset(), genotype), GT_FIELD_SEPARATOR_CHAR_);

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
        std::vector<std::string> ad_vector = Utility::char_tokenizer(genotype_parser.getFormatString(recordParser.ADOffset(), genotype), AD_FIELD_SEPARATOR_CHAR_);

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

        if (valid_record and A_allele > 0) {

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
            VariantEvidence evidence(vcf_record_count, info_evidence_opt, format_data_ptr, allele_index, allele_count);

            if (not createAddVariant(genome_name,
                                     recordParser.contigPtr(),
                                     recordParser.offset(),
                                     recordParser.reference(),
                                     allele,
                                     evidence)) {

              ExecEnv::log().error("Parsing Pf3k VCF, Problem parsing A allele CIGAR items");

            }

            ++variant_count_;

          }

        }

        if (valid_record and B_allele > 0) {

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
            VariantEvidence evidence(vcf_record_count, info_evidence_opt, format_data_ptr, allele_index, allele_count);

            if (not createAddVariant(genome_name,
                                     recordParser.contigPtr(),
                                     recordParser.offset(),
                                     recordParser.reference(),
                                     allele,
                                     evidence)) {

              ExecEnv::log().error("Parsing Pf3k VCF, Problem parsing B allele CIGAR items");

            }

            ++variant_count_;

          }

        }

      }

    }

    // Next genome name.
    ++genotype_count;

  }

  if (vcf_record_count % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("Processed :{} records, total variants: {}", vcf_record_count, variant_count_);
    ExecEnv::log().info("Contig: {}, offset: {}", recordParser.contigPtr()->contigId(), recordParser.offset());

  }

}


// Set up the genomes/contigs first rather than on-the-fly.
// Some genomes may have no variants (e.g. the model/reference genome 3D7)
// and thus these genomes/contigs would not be created on-the-fly.
void kgl::Pf3kVCFImpl::setupPopulationStructure(const std::shared_ptr<const GenomeReference> genome_db_ptr) {

  ExecEnv::log().info("setupPopulationStructure; Creating a population of {} genomes and {} contigs",
                      getGenomeNames().size(), genome_db_ptr->getMap().size());

  for (auto genome_id : getGenomeNames())  {

    std::optional<std::shared_ptr<UnphasedGenome>> genome_opt = unphased_population_ptr_->getCreateGenome(genome_id);
    if (not genome_opt) {
      // Terminate runtime.
      ExecEnv::log().critical("Could not create genome: {} in the unphased population", genome_id);

    }

    for (auto contig : genome_db_ptr->getMap()) {

      if (not genome_opt.value()->getCreateContig(contig.first)) {
        // Terminate runtime.
        ExecEnv::log().critical("Could not create contig: {} in genome: {} in the unphased population", contig.first, genome_id);

      }

    }

  }

}

bool kgl::Pf3kVCFImpl::createAddVariant(const std::string& genome_name,
                                         const std::shared_ptr<const ContigReference> contig_ptr,
                                         ContigOffset_t contig_offset,
                                         const std::string& reference_text,
                                         const std::string& alternate_text,
                                         const VariantEvidence& evidence)  {

  StringDNA5 reference_str(reference_text);
  StringDNA5 alternate_str(alternate_text);

  std::unique_ptr<const Variant> variant_ptr(std::make_unique<Variant>( contig_ptr->contigId(),
                                                                        VariantSequence::UNPHASED,
                                                                        contig_offset,
                                                                        evidence,
                                                                        std::move(reference_str),
                                                                        std::move(alternate_str)));

  return addThreadSafeVariant(std::move(variant_ptr), genome_name);

}


bool kgl::Pf3kVCFImpl::addThreadSafeVariant(std::unique_ptr<const Variant>&& variant_ptr, GenomeId_t genome) const {

  std::vector<GenomeId_t> genome_vector;
  genome_vector.emplace_back(std::move(genome));

  return unphased_population_ptr_->addVariant(std::move(variant_ptr), genome_vector);

}
