//
// Created by kellerberrin on 22/01/18.
//
#include "kgl_variant_factory_gatkvcf_impl.h"
#include "kgl_variant_factory_compound.h"


std::shared_ptr<kgl::GenomeVariant>
kgl::VcfFactory::GATKVCFImpl::readParseGATKVcfFile(const std::string& genome_name,
                                                   std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                   const std::string& vcf_file_name,
                                                   Phred_t variant_quality) {

  std::shared_ptr<GenomeVariant> genome_single_variants = kgl::GenomeVariant::emptyGenomeVariant(genome_name, genome_db_ptr);


  // Open input file.
  seqan::VcfFileIn vcfIn(seqan::toCString(vcf_file_name));

  // Attach VcfFileOut to string stream to dump record information.
  std::stringstream ss;
  seqan::VcfFileOut vcfOut(vcfIn);
  seqan::open(vcfOut, ss, seqan::Vcf());

  // Copy over header.
  seqan::VcfHeader header;
  readHeader(header, vcfIn);
  writeHeader(vcfOut, header);

  // Investigate header.
  ActiveContigMap active_contig_map;
  if (not parseVcfHeader(genome_db_ptr, header, active_contig_map, false)) {

    ExecEnv::log().error("Problem parsing header information in VCF file: {}. No variants processed.", vcf_file_name);
    ss.str("");
    writeHeader(vcfOut, header);
    ExecEnv::log().error("The VCF file header:\n{}", ss.str());
    return genome_single_variants;

  }

  // Process records.
  // Copy the file record by record.
  vcf_record_count_ = 0;
  vcf_record_error_ = 0;
  vcf_record_ignored_ = 0;
  vcf_record_rejected_ = 0;
  vcf_variant_count_ = 0;


  seqan::VcfRecord record;

  while (!seqan::atEnd(vcfIn))
  {

    size_t record_variants = 0;

    readRecord(record, vcfIn);

    ++vcf_record_count_;

    ContigId_t contig_id = seqan::toCString(contigNames(context(vcfIn))[record.rID]);
    std::shared_ptr<const ContigFeatures> contig_ptr;

    if (genome_db_ptr->getContigSequence(contig_id, contig_ptr)) {

      bool record_quality_ok;

      if (not parseVcfRecord(genome_name,
                             record,
                             contig_ptr,
                             genome_single_variants,
                             variant_quality,
                             record_quality_ok,
                             record_variants)) {

        ++vcf_record_error_;
        ss.str("");
        writeRecord(vcfOut, record);
        ExecEnv::log().error("Error parsing VCF record:\n{}", ss.str());

      }

      if (not record_quality_ok) {

        ++vcf_record_rejected_;

      }

    } else {

      ++vcf_record_ignored_;

    }


    for (size_t idx = 0; idx < record_variants; ++idx) {

      ++vcf_variant_count_;

      if (vcf_variant_count_ % VARIANT_REPORT_INTERVAL_ == 0) {

        ExecEnv::log().info("VCF file, generated: {} variants", vcf_variant_count_);

      }

    }

  }

  ExecEnv::log().info("VCF file Records; Read: {}, Rejected: {} (quality={}), Ignored: {} (no matching contig), Error: {}",
                      vcf_record_count_, vcf_record_rejected_, variant_quality, vcf_record_ignored_, vcf_record_error_);

  ExecEnv::log().info("VCF file Variants; Total generated: {}, Variant database contains :{}, Identical variants ignored: {}",
                      vcf_variant_count_, genome_single_variants->size(), vcf_variant_count_ - genome_single_variants->size());


  return genome_single_variants;

}


bool kgl::VcfFactory::GATKVCFImpl::parseVcfRecord(const std::string& genome_name,
                                                       const seqan::VcfRecord& record,
                                                       std::shared_ptr<const ContigFeatures> contig_ptr,
                                                       std::shared_ptr<GenomeVariant> genome_variants,
                                                       Phred_t variant_quality,
                                                       bool& quality_ok,
                                                       size_t& record_variants) const {


  Phred_t quality = record.qual;
  if (quality >= variant_quality) {

    if (record.beginPos < 0) {

      ExecEnv::log().error("Negative offset :{}", record.beginPos);
      return false;

    }

    quality_ok = true;
    ContigOffset_t contig_offset = static_cast<ContigOffset_t >(record.beginPos);
    std::string reference = seqan::toCString(record.ref);
    std::string alternate = seqan::toCString(record.alt);
    std::string info = toCString(record.info);
    bool result = true;

    // check sizes. both size of 1 is a SNP
    if (reference.size() == alternate.size() and alternate.size() == 1) {

      result = parseSNP(genome_name,
                        contig_ptr,
                        genome_variants,
                        quality,
                        info,
                        reference,
                        alternate,
                        contig_offset,
                        record_variants);

    } else if (reference.size() == 1 and alternate.size() > reference.size()) {

      result = parseInsert(genome_name,
                           contig_ptr,
                           genome_variants,
                           quality,
                           info,
                           reference,
                           alternate,
                           contig_offset,
                           record_variants);

    }
    else if (alternate.size() == 1 and reference.size() > alternate.size()) {

      result = parseDelete(genome_name,
                           contig_ptr,
                           genome_variants,
                           quality,
                           info,
                           reference,
                           alternate,
                           contig_offset,
                           record_variants);

    }
    else {

      ExecEnv::log().info("GATK VCF unknown variant type; reference: {}, alternate: {}", reference, alternate);
      quality_ok = false;

    }

    if (not result) {

      ExecEnv::log().error("VCF file, problem parsing GATK record.");
      return false;

    }

  } else {

    quality_ok = false;

  }

  return true;

}


bool kgl::VcfFactory::GATKVCFImpl::parseSNP(const std::string& variant_source,
                                            std::shared_ptr<const ContigFeatures> contig_ptr,
                                            std::shared_ptr<GenomeVariant> genome_variants,
                                            Phred_t quality,
                                            const std::string&, // info,
                                            const std::string& reference,
                                            const std::string& alternate,
                                            ContigOffset_t contig_offset,
                                            size_t& variant_count) const {

  // Check alternate and reference sizes.
  if (not (reference.size() == alternate.size() and reference.size() == 1)) {

    ExecEnv::log().info("GATK VCF invalid alternate and reference size for SNP; reference: {}, alternate: {}", reference, alternate);
    return false;

  }

  if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[0]) {

    ExecEnv::log().error("VCF record reference = base: {} does not match contig: {}[{}] = base: {}",
                         reference[0], contig_ptr->contigId(),
                         contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
    return false;

  }

//  std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));
  std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>("", quality));

  std::shared_ptr<SNPVariant> snp_variant_ptr(std::make_shared<SNPVariant>(variant_source,
                                                                           contig_ptr,
                                                                           contig_offset,
                                                                           quality,
                                                                           evidence_ptr,
                                                                           DNA5::convertChar(reference[0]),
                                                                           DNA5::convertChar(alternate[0])));

  variant_count += VariantFactory::addGenomeVariant(genome_variants, snp_variant_ptr); // Annotate with genome information

  return true;

}


bool kgl::VcfFactory::GATKVCFImpl::parseInsert(const std::string& variant_source,
                                               std::shared_ptr<const ContigFeatures> contig_ptr,
                                               std::shared_ptr<GenomeVariant> genome_variants,
                                               Phred_t quality,
                                               const std::string&, // info,
                                               const std::string& reference,
                                               const std::string& alternate,
                                               ContigOffset_t contig_offset,
                                               size_t& variant_count) const {

  CompoundVariantMap compound_variant_map;

  // Check alternate and reference sizes.
  if (not (reference.size() == 1 and alternate.size() > reference.size())) {

    ExecEnv::log().info("GATK VCF invalid alternate and reference size for Insert; reference: {}, alternate: {}", reference, alternate);
    return false;

  }

  // Check that the reference is valid.
  if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[0]) {

    ExecEnv::log().error("VCF record reference = base: {} does not match contig: {}[{}] = base: {}",
                         reference[0], contig_ptr->contigId(),
                         contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
    return false;

  }

  // Insert after the reference;
  ContigOffset_t insert_offset = contig_offset;
  ++insert_offset;
  ContigOffset_t alternate_offset = 1;
  ContigSize_t insert_size = alternate.size() - reference.size();
  for (size_t idx = 0; idx < insert_size; ++idx) {

//    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));
    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>("", quality));

    std::shared_ptr<InsertVariant> insert_variant_ptr(std::make_shared<InsertVariant>(variant_source,
                                                                                      contig_ptr,
                                                                                      insert_offset,
                                                                                      quality,
                                                                                      evidence_ptr,
                                                                                      contig_ptr->sequence().at(insert_offset),
                                                                                      DNA5::convertChar(alternate[alternate_offset])));


    std::pair<ContigOffset_t, std::shared_ptr<const SingleVariant>> insert_pair(insert_offset, insert_variant_ptr);
    auto result = compound_variant_map.insert(insert_pair);

    if (not result.second) {

      ExecEnv::log().error("parseInsert(), unable to insert variant with duplicate offset: {}", contig_offset);

    }

    ++alternate_offset;
    ++insert_offset;

  }

  if (compound_variant_map.size() > 1) {

    variant_count += VariantFactory::addGenomeVariant(genome_variants,CompoundInsertFactory().createCompoundVariant(compound_variant_map));

  } else if (compound_variant_map.size() == 1) {

    variant_count += VariantFactory::addGenomeVariant(genome_variants, compound_variant_map.begin()->second); // Annotate with genome information

  } else {

    ExecEnv::log().warn("parseInsert(), zero insert variants generated");

  }


  return true;

}


bool kgl::VcfFactory::GATKVCFImpl::parseDelete(const std::string& variant_source,
                                               std::shared_ptr<const ContigFeatures> contig_ptr,
                                               std::shared_ptr<GenomeVariant> genome_variants,
                                               Phred_t quality,
                                               const std::string&, // info,
                                               const std::string& reference,
                                               const std::string& alternate,
                                               ContigOffset_t contig_offset,
                                               size_t& variant_count) const {

  CompoundVariantMap compound_variant_map;

  // Check alternate and reference sizes.
  if (not (alternate.size() == 1 and reference.size() > alternate.size())) {

    ExecEnv::log().info("GATK VCF invalid alternate and reference size for Delete; reference: {}, alternate: {}", reference, alternate);
    return false;

  }

  // Check that the reference is valid.
  if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[0]) {

    ExecEnv::log().error("VCF record reference = base: {} does not match contig: {}[{}] = base: {}",
                         reference[0], contig_ptr->contigId(),
                         contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
    return false;

  }

  // Delete after the reference;
  ContigOffset_t delete_offset = contig_offset;
  ++delete_offset;
  ContigOffset_t reference_offset = 1;
  ContigSize_t delete_size = reference.size() - alternate.size();

  for (size_t idx = 0; idx < delete_size; ++idx) {

    if (DNA5::convertToChar(contig_ptr->sequence().at(delete_offset)) != reference[reference_offset]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match contig: {}[{}] = base: {}",
                           reference_offset, reference[reference_offset], contig_ptr->contigId(),
                           contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(delete_offset)));
      return false;

    }

//    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));
    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>("", quality));

    std::shared_ptr<DeleteVariant> delete_variant_ptr(std::make_shared<DeleteVariant>(variant_source,
                                                                                      contig_ptr,
                                                                                      delete_offset,
                                                                                      quality,
                                                                                      evidence_ptr,
                                                                                      contig_ptr->sequence().at(delete_offset)));

    std::pair<ContigOffset_t, std::shared_ptr<const SingleVariant>> insert_pair(delete_offset, delete_variant_ptr);
    auto result = compound_variant_map.insert(insert_pair);

    if (not result.second) {

      ExecEnv::log().error("parseDelete(), unable to delete variant with duplicate offset: {}", delete_offset);

    }

    ++reference_offset;
    ++delete_offset;

  }

  if (compound_variant_map.size() > 1) {

    variant_count += VariantFactory::addGenomeVariant(genome_variants,CompoundDeleteFactory().createCompoundVariant(compound_variant_map));

  } else if (compound_variant_map.size() == 1) {

    variant_count += VariantFactory::addGenomeVariant(genome_variants, compound_variant_map.begin()->second); // Annotate with genome information

  } else {

    ExecEnv::log().error("parseDelete(), zero delete variants generated");

  }


  return true;

}
