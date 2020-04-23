//
// Created by kellerberrin on 25/02/18.
//

#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_pf3k_impl.h"
#include "kgl_variant_vcf.h"


namespace kgl = kellerberrin::genome;


// Process VCF header information.
void kgl::Pf3kVCFImpl::processVCFHeader(const VcfHeaderInfo& header_info) {

  // Investigate header.
  ActiveContigMap active_contig_map;
  if (not ParseVCFMiscImpl::parseVcfHeader(genome_db_ptr_, header_info, active_contig_map, false)) {

    ExecEnv::log().error("Problem parsing header information in VCF file. No variants processed.");

  }

  // Pre-initializes the UnphasedPopulation with a list of genomes/contigs.
  // Since some genomes may not have a variant (3D7).
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

  // Parse the info fields into an array and assign to a shared ptr.
  VCFInfoField info_key_value_map(record.info);  // Each vcf record.

  // For performance reasons the info field is moved here - don't it reference again.
  auto mutable_info = const_cast<std::string&>(record.info);
  std::shared_ptr<std::string> info_ptr = std::make_shared<std::string>(std::move(mutable_info));

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
            std::shared_ptr<VariantEvidence> evidence_ptr(std::make_shared<CountEvidence>(info_ptr,
                                                                                          ref_count,
                                                                                          alt_count,
                                                                                          DP_value,
                                                                                          GQ_value,
                                                                                          recordParser.quality(),
                                                                                          vcf_record_count));

            if (not createAddVariant(genome_name,
                                     recordParser.contigPtr(),
                                     recordParser.offset(),
                                     recordParser.reference(),
                                     allele,
                                     evidence_ptr)) {

              ExecEnv::log().error("Parsing Pf3k VCF, Problem parsing A allele CIGAR items");

            }

            ++variant_count_;

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
            std::shared_ptr<VariantEvidence> evidence_ptr(std::make_shared<CountEvidence>(info_ptr,
                                                                                          ref_count,
                                                                                          alt_count,
                                                                                          DP_value,
                                                                                          GQ_value,
                                                                                          recordParser.quality(),
                                                                                          vcf_record_count));


            if (not createAddVariant(genome_name,
                                     recordParser.contigPtr(),
                                     recordParser.offset(),
                                     recordParser.reference(),
                                     allele,
                                     evidence_ptr)) {

              ExecEnv::log().error("Parsing Pf3k VCF, Problem parsing B allele CIGAR items");

            }

            ++variant_count_;

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

bool kgl::Pf3kVCFImpl::addThreadSafeGenomeVariant(const std::shared_ptr<const Variant>& variant_ptr) {

  std::scoped_lock<std::mutex> auto_mutex(add_variant_mutex_); // Write Locked

  std::shared_ptr<UnphasedGenome> genome;

  if (not unphased_population_ptr_->getCreateGenome(variant_ptr->genomeId(), genome)) {

    ExecEnv::log().error("ParseVCFImpl::addThreadSafeGenomeVariant; Could not add/create genome: {}", variant_ptr->genomeId());
    return false;

  }

  if (not genome->addVariant(variant_ptr)) { // thread safe

    ExecEnv::log().error("ParseVCFImpl::addThreadSafeGenomeVariant; Could not add variant to genome: {}", variant_ptr->genomeId());
    return false;

  }

  return true;

}


// Set up the genomes/contigs first rather than on-the-fly.
// Some genomes may have no variants (e.g. the model/reference genome 3D7)
// and thus these genomes/contigs would not be created on-the-fly.
void kgl::Pf3kVCFImpl::setupPopulationStructure(const std::shared_ptr<const GenomeDatabase> genome_db_ptr) {

  std::scoped_lock<std::mutex> lock(add_variant_mutex_);

  ExecEnv::log().info("setupPopulationStructure; Creating a population of {} genomes and {} contigs",
                      getGenomeNames().size(), genome_db_ptr->getMap().size());

  for (auto genome_id : getGenomeNames())  {

    std::shared_ptr<UnphasedGenome> genome_ptr = nullptr;

    if (not unphased_population_ptr_->getCreateGenome(genome_id, genome_ptr)) {

      ExecEnv::log().critical("Could not create genome: {} in the unphased population", genome_id);

    }

    for (auto contig : genome_db_ptr->getMap()) {

      std::shared_ptr<UnphasedContig> contig_ptr = nullptr;

      if (not genome_ptr->getCreateContig(contig.first, contig_ptr)) {

        ExecEnv::log().critical("Could not create contig: {} in genome: {} in the unphased population", contig.first, genome_id);

      }

    }

  }

}

bool kgl::Pf3kVCFImpl::createAddVariant(const std::string& genome_name,
                                         const std::shared_ptr<const ContigFeatures> contig_ptr,
                                         ContigOffset_t contig_offset,
                                         const std::string& reference_text,
                                         const std::string& alternate_text,
                                         const std::shared_ptr<const VariantEvidence> evidence_ptr)  {

  StringDNA5 reference_str(reference_text);
  StringDNA5 alternate_str(alternate_text);

  std::shared_ptr<const Variant> variant_ptr(std::make_shared<VCFVariant>(genome_name,
                                                                          contig_ptr->contigId(),
                                                                          VariantSequence::UNPHASED,
                                                                          contig_offset,
                                                                          evidence_ptr,
                                                                          std::move(reference_str),
                                                                          std::move(alternate_str)));

  return unphased_population_ptr_->addThreadSafeVariant(variant_ptr);

}
