//
// Created by kellerberrin on 22/01/18.
//


#include "kgl_variant_factory_fbvcf_impl.h"
#include "kgl_variant_factory_compound.h"

namespace kgl = kellerberrin::genome;
namespace bt = boost;


std::shared_ptr<kgl::GenomeVariant>
kgl::VcfFactory::FreeBayesVCFImpl::readParseFreeBayesVcfFile(const std::string& genome_name,
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
  if (not parseVcfHeader(genome_db_ptr, header, active_contig_map, true)) {

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

bool kgl::VcfFactory::FreeBayesVCFImpl::parseVcfRecord(const std::string& genome_name,
                                                  const seqan::VcfRecord& record,
                                                  std::shared_ptr<const ContigFeatures> contig_ptr,
                                                  std::shared_ptr<GenomeVariant> genome_variants,
                                                  Phred_t variant_quality,
                                                  bool& quality_ok,
                                                  size_t& record_variants) const {

  std::string info = toCString(record.info);
  // assumes input "key_1=value_1; ...;key_n=value_n"
  std::map<std::string, std::string> info_key_value_map;
  if (not tokenizeVcfInfoKeyValues(info, info_key_value_map)) {

    ExecEnv::log().error("Unable to parse VCF INFO: {}", info);
    return false;

  }

  std::vector<std::pair<char, size_t>> parsed_cigar;
  size_t reference_size;
  size_t alternate_size;
  std::string cigar;
  Phred_t quality = record.qual;

  auto result_cigar = info_key_value_map.find(ID_CIGAR_VALUE_);

  if (result_cigar != info_key_value_map.end()) {

    cigar = result_cigar->second;
    if (not parseCigar(cigar, reference_size, alternate_size, parsed_cigar)) {

      return false;

    }

  } else {

    ExecEnv::log().error("VCF factory; cigar field: {} not found in info: {}", ID_CIGAR_VALUE_, info);
    ExecEnv::log().error("VCF file should conform to 'freebayes' format");
    return false;

  }

  if (quality >= variant_quality) {

    quality_ok = true;
    std::string reference = seqan::toCString(record.ref);
    std::string alternate = seqan::toCString(record.alt);
    if (record.beginPos < 0) {

      ExecEnv::log().error("");

    }

    // check sizes.
    if (reference.size() != reference_size) {

      ExecEnv::log().error("VCF factory; reference: {} size: {} does not match cigar: {} size: {}",
                           reference, reference.size(), cigar, reference_size);
      return false;

    }

    if (alternate.size() != alternate_size) {

      ExecEnv::log().error("VCF factory; alternative: {} size: {} does not match cigar: {} size: {}",
                           alternate, alternate.size(), cigar, alternate_size);
      return false;

    }

    ContigOffset_t contig_offset = static_cast<ContigOffset_t >(record.beginPos);
    size_t reference_index = 0;
    size_t alternate_index = 0;
    bool result = false;

    // Generate variants.
    for (auto cigar_item : parsed_cigar) {

      switch(cigar_item.first) {


        case 'M':
          result = parseCheck(cigar_item.second,
                              contig_ptr,
                              reference,
                              alternate,
                              reference_index,
                              alternate_index,
                              contig_offset);
          break;

        case 'X':
          result = parseSNP(cigar_item.second,
                            genome_name,
                            contig_ptr,
                            genome_variants,
                            quality,
                            info,
                            reference,
                            alternate,
                            reference_index,
                            alternate_index,
                            contig_offset,
                            record_variants);
          break;

        case 'I':
          result = parseInsert(cigar_item.second,
                               genome_name,
                               contig_ptr,
                               genome_variants,
                               quality,
                               info,
                               alternate,
                               contig_offset,
                               alternate_index,
                               record_variants);
          break;

        case 'D':
          result = parseDelete(cigar_item.second,
                               genome_name,
                               contig_ptr,
                               genome_variants,
                               quality,
                               info,
                               reference,
                               reference_index,
                               contig_offset,
                               record_variants);
          break;

      }

      if (not result) {

        ExecEnv::log().error("VCF file, problem parsing cigar element {}:{}", cigar_item.second, cigar_item.first);
        return false;

      }

    }

  } else {

    quality_ok = false;

  }

  return true;

}


bool kgl::VcfFactory::FreeBayesVCFImpl::parseCheck(size_t cigar_count,
                                              std::shared_ptr<const ContigFeatures> contig_ptr,
                                              const std::string& reference,
                                              const std::string& alternate,
                                              size_t& reference_index,
                                              size_t& alternate_index,
                                              ContigOffset_t& contig_offset) const {

  for (size_t idx = 0; idx < cigar_count; ++idx) {

    if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[reference_index]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match contig: {}[{}] = base: {}",
                           reference_index, reference[reference_index], contig_ptr->contigId(),
                           contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
      return false;

    }

    if (alternate[alternate_index] != reference[reference_index]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match VCF record alternate[{}] = base: {}",
                           reference_index, reference[reference_index], alternate_index, alternate[alternate_index]);
      return false;

    }

    ++reference_index;
    ++alternate_index;
    ++contig_offset;

  }

  return true;

}


bool kgl::VcfFactory::FreeBayesVCFImpl::parseSNP(size_t cigar_count,
                                            const std::string& variant_source,
                                            std::shared_ptr<const ContigFeatures> contig_ptr,
                                            std::shared_ptr<GenomeVariant> genome_variants,
                                            Phred_t quality,
                                            const std::string&, // info,
                                            const std::string& reference,
                                            const std::string& alternate,
                                            size_t& reference_index,
                                            size_t& alternate_index,
                                            ContigOffset_t& contig_offset,
                                            size_t& variant_count) const {

  for (size_t idx = 0; idx < cigar_count; ++idx) {

    if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[reference_index]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match contig: {}[{}] = base: {}",
                           reference_index, reference[reference_index], contig_ptr->contigId(),
                           contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
      return false;

    }

//    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));
    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>("", quality));

    std::shared_ptr<SNPVariant> snp_variant_ptr(std::make_shared<SNPVariant>(variant_source,
                                                                             contig_ptr,
                                                                             contig_offset,
                                                                             quality,
                                                                             evidence_ptr,
                                                                             DNA5::convertChar(reference[reference_index]),
                                                                             DNA5::convertChar(alternate[alternate_index])));

    variant_count += VariantFactory::addGenomeVariant(genome_variants, snp_variant_ptr); // Annotate with genome information

    ++reference_index;
    ++alternate_index;
    ++contig_offset;

  }

  return true;

}


bool kgl::VcfFactory::FreeBayesVCFImpl::parseInsert(size_t cigar_count,
                                               const std::string& variant_source,
                                               std::shared_ptr<const ContigFeatures> contig_ptr,
                                               std::shared_ptr<GenomeVariant> genome_variants,
                                               Phred_t quality,
                                               const std::string&, // info,
                                               const std::string& alternate,
                                               ContigOffset_t contig_offset,
                                               size_t& alternate_index,
                                               size_t& variant_count) const {

  CompoundVariantMap compound_variant_map;

  for (size_t idx = 0; idx < cigar_count; ++idx) {

//    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));
    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>("", quality));

    std::shared_ptr<InsertVariant> insert_variant_ptr(std::make_shared<InsertVariant>(variant_source,
                                                                                      contig_ptr,
                                                                                      contig_offset,
                                                                                      quality,
                                                                                      evidence_ptr,
                                                                                      contig_ptr->sequence().at(contig_offset),
                                                                                      DNA5::convertChar(alternate[alternate_index])));


    std::pair<ContigOffset_t, std::shared_ptr<const SingleVariant>> insert_pair(contig_offset, insert_variant_ptr);
    auto result = compound_variant_map.insert(insert_pair);

    if (not result.second) {

      ExecEnv::log().error("parseInsert(), unable to insert variant with duplicate offset: {}", contig_offset);

    }

    ++alternate_index;
    ++contig_offset;

  }

  if (compound_variant_map.size() > 1) {

    variant_count += VariantFactory::addGenomeVariant(genome_variants,CompoundInsertFactory().createCompoundVariant(compound_variant_map));

  } else if (compound_variant_map.size() == 1) {

    variant_count += VariantFactory::addGenomeVariant(genome_variants, compound_variant_map.begin()->second); // Annotate with genome information

  } else {

    ExecEnv::log().error("parseInsert(), cigar size: {} but zero insert variant generated", cigar_count);

  }


  return true;

}

bool kgl::VcfFactory::FreeBayesVCFImpl::parseDelete(size_t cigar_count,
                                               const std::string& variant_source,
                                               std::shared_ptr<const ContigFeatures> contig_ptr,
                                               std::shared_ptr<GenomeVariant> genome_variants,
                                               Phred_t quality,
                                               const std::string&, // info,
                                               const std::string& reference,
                                               size_t& reference_index,
                                               ContigOffset_t& contig_offset,
                                               size_t& variant_count) const {

  CompoundVariantMap compound_variant_map;

  for (size_t idx = 0; idx < cigar_count; ++idx) {

    if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[reference_index]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match contig: {}[{}] = base: {}",
                           reference_index, reference[reference_index], contig_ptr->contigId(),
                           contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
      return false;

    }

//    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));
    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>("", quality));


    std::shared_ptr<DeleteVariant> delete_variant_ptr(std::make_shared<DeleteVariant>(variant_source,
                                                                                      contig_ptr,
                                                                                      contig_offset,
                                                                                      quality,
                                                                                      evidence_ptr,
                                                                                      contig_ptr->sequence().at(contig_offset)));

    std::pair<ContigOffset_t, std::shared_ptr<const SingleVariant>> insert_pair(contig_offset, delete_variant_ptr);
    auto result = compound_variant_map.insert(insert_pair);

    if (not result.second) {

      ExecEnv::log().error("parseDelete(), unable to insert variant with duplicate offset: {}", contig_offset);

    }

    ++reference_index;
    ++contig_offset;

  }

  if (compound_variant_map.size() > 1) {

    variant_count += VariantFactory::addGenomeVariant(genome_variants,CompoundDeleteFactory().createCompoundVariant(compound_variant_map));

  } else if (compound_variant_map.size() == 1) {

    variant_count += VariantFactory::addGenomeVariant(genome_variants, compound_variant_map.begin()->second); // Annotate with genome information

  } else {

    ExecEnv::log().error("parseDelete(), cigar size: {} but zero insert variant generated", cigar_count);

  }


  return true;

}

