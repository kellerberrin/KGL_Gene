//
// Created by kellerberrin on 26/12/17.
//


#include "kgl_utility.h"
#include "kgl_sam_process.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_single.h"
#include "kgl_variant_factory_single.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <seqan/vcf_io.h>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory::VcfFileImpl does all the heavy lifting using 3rd a party library. In this case; Seqan.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ActiveContigMap = std::map<kgl::ContigId_t, kgl::ContigSize_t>;

class kgl::VcfFactory::VcfFileImpl {

public:

  VcfFileImpl() = default;
  ~VcfFileImpl() = default;

  std::shared_ptr<GenomeVariant> readParseVcfFile(const std::string& genome_name,
                                             std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                             const std::string& vcf_file_name,
                                             Phred_t variant_quality);

private:

  constexpr static const char* HEADER_CONTIG_KEY_ = "CONTIG";
  constexpr static const char* ID_KEY_ = "ID";
  constexpr static const char* CONTIG_LENGTH_KEY_ = "LENGTH";
  constexpr static const char* HEADER_INFO_KEY_ = "INFO";
  constexpr static const char* ID_CIGAR_VALUE_ = "CIGAR";
  constexpr static const char* ID_READ_DEPTH_ = "DPB";
  constexpr static const char* ID_PROPORTION_ = "AF";

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 5000;

  size_t vcf_record_count_;
  size_t vcf_record_ignored_;
  size_t vcf_record_error_;
  size_t vcf_record_rejected_;

  bool parseVcfHeader(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                      const seqan::VcfHeader& header,
                      ActiveContigMap& active_contig_map) const;

  bool parseVcfRecord(const std::string& genome_name,
                      const seqan::VcfRecord& record,
                      ActiveContigMap& active_contig_map,
                      std::shared_ptr<const ContigFeatures> contig_ptr,
                      std::shared_ptr<GenomeVariant> genome_variants,
                      Phred_t variant_quality,
                      bool& quality_ok) const;

  // Parse 1M ... XM in the cigar.
  bool parseCheck(size_t cigar_count,
                  std::shared_ptr<const ContigFeatures> contig_ptr,
                  const std::string& reference,
                  const std::string& alternate,
                  size_t& reference_index,
                  size_t& alternate_index,
                  ContigOffset_t& contig_offset) const;

  // Parse 1X ... XX in the cigar.

  bool parseSNP(size_t cigar_count,
                const std::string& variant_source,
                std::shared_ptr<const ContigFeatures> contig_ptr,
                std::shared_ptr<GenomeVariant> genome_variants,
                Phred_t quality,
                const std::string& info,
                const std::string& reference,
                const std::string& alternate,
                size_t& reference_index,
                size_t& alternate_index,
                ContigOffset_t& contig_offset,
                size_t& variant_count) const;

  // Parse 1I ... XI in the cigar.
  bool parseInsert(size_t cigar_count,
                   const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& alternate,
                   size_t& alternate_index,
                   ContigOffset_t& contig_offset,
                   size_t& variant_count) const;

  // Parse 1D ... XD in the cigar.
  bool parseDelete(size_t cigar_count,
                   const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& reference,
                   size_t& reference_index,
                   ContigOffset_t& contig_offset,
                   size_t& variant_count) const;

  // assumes input "<key_1=value_1, ...,key_n=value_n>"
  bool tokenizeVcfHeaderKeyValues(const std::string& key_value_text,
                                  std::map<std::string, std::string>& key_value_pairs) const;

// assumes input "key_1=value_1; ...;key_n=value_n"
  bool tokenizeVcfInfoKeyValues(const std::string& key_value_text,
                                std::map<std::string, std::string>& key_value_map) const;

  bool parseCigar(const std::string& cigar,
                  size_t& check_reference_size,
                  size_t& check_alternate_size,
                  std::vector<std::pair<char, size_t>>& parsed_cigar) const;

};


std::shared_ptr<kgl::GenomeVariant>
kgl::VcfFactory::VcfFileImpl::readParseVcfFile(const std::string& genome_name,
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
  if (not parseVcfHeader(genome_db_ptr, header, active_contig_map)) {

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

  seqan::VcfRecord record;

  while (!seqan::atEnd(vcfIn))
  {

    readRecord(record, vcfIn);

    ++vcf_record_count_;

    ContigId_t contig_id = seqan::toCString(contigNames(context(vcfIn))[record.rID]);
    std::shared_ptr<const ContigFeatures> contig_ptr;
    if (genome_db_ptr->getContigSequence(contig_id, contig_ptr)) {

      bool record_quality_ok;
      if (not parseVcfRecord(genome_name,
                             record,
                             active_contig_map,
                             contig_ptr,
                             genome_single_variants,
                             variant_quality,
                             record_quality_ok)) {

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

  }

  ExecEnv::log().info("VCF file; Read: {}, Rejected: {} (quality={}), Ignored: {} (no matching contig), Error: {} records.",
                      vcf_record_count_, vcf_record_rejected_, variant_quality, vcf_record_ignored_, vcf_record_error_);
  ExecEnv::log().info("VCF file: {} generated: {} raw variants.", vcf_file_name, genome_single_variants->size());

  return genome_single_variants;

}

bool kgl::VcfFactory::VcfFileImpl::parseVcfRecord(const std::string& genome_name,
                                                  const seqan::VcfRecord& record,
                                                  ActiveContigMap& active_contig_map,
                                                  std::shared_ptr<const ContigFeatures> contig_ptr,
                                                  std::shared_ptr<GenomeVariant> genome_variants,
                                                  Phred_t variant_quality,
                                                  bool& quality_ok) const {

  /*
  record.info
  record.id
  record.alt
  record.beginPos
  record.filter
  record.format
  record.genotypeInfos
  record.INVALID_POS
  record.INVALID_REFID
  record.qual
  record.ref
  record.rID
  record.MISSING_QUAL() */

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
    size_t variant_total = 0;
    bool result = false;


    // Generate variants.
    for (auto cigar_item : parsed_cigar) {

      size_t variant_count = 0;

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
                            variant_count);
          break;

        case 'I':
          result = parseInsert(cigar_item.second,
                               genome_name,
                               contig_ptr,
                               genome_variants,
                               quality,
                               info,
                               alternate,
                               alternate_index,
                               contig_offset,
                               variant_count);
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
                               variant_count);
          break;

      }

      if (not result) {

        ExecEnv::log().error("VCF file, problem parsing cigar element {}:{}", cigar_item.second, cigar_item.first);
        return false;

      }

      for (size_t idx = 0; idx < variant_count; ++idx) {

        ++variant_total;

        if (variant_total % VARIANT_REPORT_INTERVAL_ == 0) {

          ExecEnv::log().info("VCF file, generated: {} variants", variant_total);

        }

      }

    }

  } else {

    quality_ok = false;

  }

  return true;

}


bool kgl::VcfFactory::VcfFileImpl::parseCheck(size_t cigar_count,
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


bool kgl::VcfFactory::VcfFileImpl::parseSNP(size_t cigar_count,
                                            const std::string& variant_source,
                                            std::shared_ptr<const ContigFeatures> contig_ptr,
                                            std::shared_ptr<GenomeVariant> genome_variants,
                                            Phred_t quality,
                                            const std::string& info,
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

    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));

    SNPVariant snp_variant(variant_source,
                           contig_ptr,
                           contig_offset,
                           quality,
                           evidence_ptr,
                           DNA5::convertChar(reference[reference_index]),
                           DNA5::convertChar(alternate[alternate_index]));

    variant_count += VariantFactory::addSingleVariant(genome_variants, snp_variant); // Annotate with genome information

    ++reference_index;
    ++alternate_index;
    ++contig_offset;

  }

  return true;

}


bool kgl::VcfFactory::VcfFileImpl::parseInsert(size_t cigar_count,
                                               const std::string& variant_source,
                                               std::shared_ptr<const ContigFeatures> contig_ptr,
                                               std::shared_ptr<GenomeVariant> genome_variants,
                                               Phred_t quality,
                                               const std::string& info,
                                               const std::string& alternate,
                                               size_t& alternate_index,
                                               ContigOffset_t& contig_offset,
                                               size_t& variant_count) const {

  for (size_t idx = 0; idx < cigar_count; ++idx) {

    ++alternate_index;

  }

  return true;

}

bool kgl::VcfFactory::VcfFileImpl::parseDelete(size_t cigar_count,
                                               const std::string& variant_source,
                                               std::shared_ptr<const ContigFeatures> contig_ptr,
                                               std::shared_ptr<GenomeVariant> genome_variants,
                                               Phred_t quality,
                                               const std::string& info,
                                               const std::string& reference,
                                               size_t& reference_index,
                                               ContigOffset_t& contig_offset,
                                               size_t& variant_count) const {


  for (size_t idx = 0; idx < cigar_count; ++idx) {

    if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[reference_index]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match contig: {}[{}] = base: {}",
                           reference_index, reference[reference_index], contig_ptr->contigId(),
                           contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
      return false;

    }

    ++reference_index;
    ++contig_offset;

  }

  return true;

}


bool kgl::VcfFactory::VcfFileImpl::parseCigar(const std::string& cigar,
                                              size_t& check_reference_size,
                                              size_t& check_alternate_size,
                                              std::vector<std::pair<char, size_t>>& parsed_cigar) const {

  parsed_cigar.clear();
  check_reference_size = 0;
  check_alternate_size = 0;
  auto it = cigar.begin();
  char cigar_code;
  std::string cigar_size;
  while(it != cigar.end()) {

    cigar_size.clear();
    while(isdigit(*it) and it != cigar.end()) {

      cigar_size += *it;
      ++it;

    }

    if (cigar_size.empty()) {

      ExecEnv::log().error("VCF factory; cigar: {} contains no cigar size", cigar);
      return false;

    }

    size_t size = std::atoll(cigar_size.c_str());

    if (it != cigar.end()) {

      switch (*it) {

        case 'I':
          check_alternate_size += size;
          cigar_code = *it;
          break;

        case 'D':
          check_reference_size += size;
          cigar_code = *it;
          break;

        case 'X':
        case 'M':
          check_alternate_size += size;
          check_reference_size += size;
          cigar_code = *it;
          break;

        default:
          ExecEnv::log().error("VCF factory; cigar: {} contains unexpected cigar code: {}", cigar, *it);
          return false;

      }

      parsed_cigar.push_back(std::pair<char, size_t>(cigar_code, size));
      ++it;

    } else {

      ExecEnv::log().error("VCF factory; cigar: {} unexpected end; no terminating cigar code.", cigar);
      return false;

    }

  }

  return true;

}


bool kgl::VcfFactory::VcfFileImpl::parseVcfHeader(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                  const seqan::VcfHeader& header,
                                                  ActiveContigMap& active_contig_map) const {

  active_contig_map.clear();
  bool has_cigar = false;

  for (size_t idx = 0; idx != seqan::length(header); ++idx) {

    std::string key = seqan::toCString(header[idx].key);
    std::transform(key.begin(), key.end(), key.begin(), ::toupper);

    if (key == HEADER_CONTIG_KEY_) {

      std::map<std::string, std::string> item_map;
      std::string value = seqan::toCString(header[idx].value);
      tokenizeVcfHeaderKeyValues(value, item_map);

      auto id_result = item_map.find(ID_KEY_);
      auto length_result = item_map.find(CONTIG_LENGTH_KEY_);

      if (id_result != item_map.end() and length_result != item_map.end()) {

        ContigId_t contig_id = id_result->second;
        ContigSize_t contig_size = std::atoll(length_result->second.c_str());

        std::shared_ptr<const ContigFeatures> contig_ptr;
        if (genome_db_ptr->getContigSequence(contig_id, contig_ptr)) {

          if (contig_ptr->sequence().length() == contig_size) {

            active_contig_map[contig_id] = contig_size;

          } else {

            ExecEnv::log().warn("VCF header. VCF contig: {} size: {} not equal to the Genome database contig size: {}",
                                contig_id, contig_size, contig_ptr->sequence().length());
            ExecEnv::log().warn("Check that the VCF fasta file is compatible with the Genome database fasta file.");
          }

        } else {

          ExecEnv::log().warn("VCF header. VCF contig: {} size: {} not found in the genome database. ",
                              contig_id, contig_size);

        }

      }

    } else if (key == HEADER_INFO_KEY_) {

      std::map<std::string, std::string> item_map;
      std::string value = seqan::toCString(header[idx].value);
      tokenizeVcfHeaderKeyValues(value, item_map);

      auto id_result = item_map.find(ID_KEY_);

      if (id_result != item_map.end()) {

        std::string value = id_result->second;
        std::transform(value.begin(), value.end(), value.begin(), ::toupper);
        if (value == ID_CIGAR_VALUE_) {

          has_cigar = true;

        }

      }

    }

  }

  if (active_contig_map.size() == genome_db_ptr->contigCount()) {

    ExecEnv::log().info("Genome database contigs: {} and VCF contigs: {} all match.", genome_db_ptr->contigCount(), active_contig_map.size());

  } else {

    ExecEnv::log().warn("VCF file contig count: {}. Genome database contig count: {}", active_contig_map.size(), genome_db_ptr->contigCount());
    ExecEnv::log().warn("No variants will be generated for the missing contigs. The missing contigs are:");

    for (auto contig : genome_db_ptr->getMap()) {

      auto result = active_contig_map.find(contig.first);

      if (result == active_contig_map.end()) {

        ExecEnv::log().warn("Contig: {} present in genome database, missing from VCF file.", contig.first);

      }

    }

  }

  if (not has_cigar) {

    ExecEnv::log().info("This VCF file does not define a 'CIGAR' field; this field is required to parse the VCF files");
    ExecEnv::log().info("See the VCF format of the variant caller 'freebayes' for more information.");
    return false;

  }

  return true;

}

// assumes input "<key_1=value_1, ...,key_n=value_n>"
bool kgl::VcfFactory::VcfFileImpl::tokenizeVcfHeaderKeyValues(const std::string& key_value_text,
                                                              std::map<std::string, std::string>& key_value_map) const {

  key_value_map.clear();
  std::string copy_key_value_text = key_value_text;
  // Remove the angle backets from the value.
  bt::erase_all(copy_key_value_text, "<");
  bt::erase_all(copy_key_value_text, ">");


  bt::tokenizer<bt::escaped_list_separator<char>> tokenize(copy_key_value_text);
  for(auto iter = tokenize.begin(); iter != tokenize.end(); ++iter) {

    std::vector<std::string> item_vec;
    bt::char_separator<char> item_key_sep("=");
    bt::tokenizer<bt::char_separator<char>> tokenize_item(*iter, item_key_sep);
    for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

      item_vec.push_back(seqan::toCString(*iter_item));

    }

    if (item_vec.size() >= 2) {

      std::string item_key = item_vec[0];
      std::transform(item_key.begin(), item_key.end(), item_key.begin(), ::toupper);
      key_value_map[item_key] = item_vec[1];

    } else {

      ExecEnv::log().warn("Problem parsing item: {} in VCF header line: {}, expected <key=value,..> pairs", *iter, key_value_text);
      return false;

    }

  }

  return true;

}

// assumes input "key_1=value_1; ...;key_n=value_n"
bool kgl::VcfFactory::VcfFileImpl::tokenizeVcfInfoKeyValues(const std::string& key_value_text,
                                                              std::map<std::string, std::string>& key_value_map) const {

  key_value_map.clear();

  bt::char_separator<char> item_key_sep(";");
  bt::tokenizer<bt::char_separator<char>> tokenize(key_value_text, item_key_sep);
  for(auto iter = tokenize.begin(); iter != tokenize.end(); ++iter) {

    std::vector<std::string> item_vec;
    bt::char_separator<char> item_key_sep("=");
    bt::tokenizer<bt::char_separator<char>> tokenize_item(*iter, item_key_sep);
    for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

      item_vec.push_back(seqan::toCString(*iter_item));

    }

    if (item_vec.size() >= 2) {

      std::string item_key = item_vec[0];
      std::transform(item_key.begin(), item_key.end(), item_key.begin(), ::toupper);
      key_value_map[item_key] = item_vec[1];

    } else {

      ExecEnv::log().warn("Problem parsing item: {} in VCF Record line: {}, expected 'key=value;...;key=value' pairs", *iter, key_value_text);
      return false;

    }

  }

  return true;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto VcfFactory::VcfFileImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::VcfFactory::VcfFactory() : vcf_file_impl_ptr_(std::make_unique<kgl::VcfFactory::VcfFileImpl>()) {}
kgl::VcfFactory::~VcfFactory() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.


std::shared_ptr<kgl::GenomeVariant> kgl::VcfFactory::readParseVcf(const std::string& genome_name,
                                                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                                  const std::string& vcf_file_name,
                                                                  Phred_t variant_quality) const {

  return vcf_file_impl_ptr_->readParseVcfFile(genome_name,
                                              genome_db_ptr,
                                              vcf_file_name,
                                              variant_quality);

}

