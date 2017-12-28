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
                                             Phred_t read_quality,
                                             NucleotideReadCount_t min_read_count,
                                             double min_proportion) const;

private:

  constexpr static const char* HEADER_CONTIG_KEY_ = "CONTIG";
  constexpr static const char* ID_KEY_ = "ID";
  constexpr static const char* CONTIG_LENGTH_KEY_ = "LENGTH";
  constexpr static const char* HEADER_INFO_KEY_ = "INFO";
  constexpr static const char* ID_CIGAR_VALUE_ = "CIGAR";
  constexpr static const char* ID_READ_DEPTH_ = "DPB";
  constexpr static const char* ID_PROPORTION_ = "AF";

  bool parseVcfHeader(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                      const seqan::VcfHeader& header,
                      ActiveContigMap& active_contig_map) const;

  bool parseVcfRecord(const std::string& genome_name,
                      const seqan::VcfRecord& record,
                      ActiveContigMap& active_contig_map,
                      std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                      std::shared_ptr<GenomeVariant> genome_variants,
                      Phred_t read_quality,
                      NucleotideReadCount_t min_read_count,
                      double min_proportion) const;

  // assumes input "<key_1=value_1, ...,key_n=value_n>"
  bool tokenizeVcfHeaderKeyValues(const std::string& key_value_text,
                                  std::map<std::string, std::string>& key_value_pairs) const;

// assumes input "key_1=value_1; ...;key_n=value_n"
  bool tokenizeVcfInfoKeyValues(const std::string& key_value_text,
                                std::map<std::string, std::string>& key_value_map) const;

  bool parseCigar(const std::string& cigar,
                  size_t& total_cigar_size,
                  std::vector<std::pair<char, size_t>>& parsed_cigar) const;

};


std::shared_ptr<kgl::GenomeVariant>
kgl::VcfFactory::VcfFileImpl::readParseVcfFile(const std::string& genome_name,
                                               std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                               const std::string& vcf_file_name,
                                               Phred_t read_quality,
                                               NucleotideReadCount_t min_read_count,
                                               double min_proportion) const {

  std::shared_ptr<GenomeVariant> genome_single_variants = kgl::GenomeVariant::emptyGenomeVariant(genome_name, genome_db_ptr);


  // Open input file.
  seqan::VcfFileIn vcfIn(seqan::toCString(vcf_file_name));

  // Attach to string stream.
  std::stringstream ss;
  seqan::VcfFileOut vcfOut(vcfIn);
  seqan::open(vcfOut, std::cout, seqan::Vcf());

  // Copy over header.
  seqan::VcfHeader header;
  readHeader(header, vcfIn);
  writeHeader(vcfOut, header);

  ActiveContigMap active_contig_map;
  if (not parseVcfHeader(genome_db_ptr, header, active_contig_map)) {

    ExecEnv::log().error("Problem parsing header information in VCF file: {}. No variants processed.", vcf_file_name);
    ss.str("");
    writeHeader(vcfOut, header);
    ExecEnv::log().error("The VCF file header:\n{}", ss.str());
    return genome_single_variants;

  }

  // Copy the file record by record.
  size_t vcf_record_count = 0;
  seqan::VcfRecord record;

  while (!seqan::atEnd(vcfIn))
  {

    readRecord(record, vcfIn);

    ++vcf_record_count;
    if (vcf_record_count < 1000000000) {

//      writeRecord(vcfOut, record);

      if (not parseVcfRecord(genome_name,
                             record,
                             active_contig_map,
                             genome_db_ptr,
                             genome_single_variants,
                             read_quality,
                             min_read_count,
                             min_proportion)) {

      }

    }

  }

  ExecEnv::log().info("Read: {} records from VCF file: {}", vcf_record_count, vcf_file_name);
  ExecEnv::log().info("VCF file {} has: {} raw variants", vcf_file_name, genome_single_variants->size());

  return genome_single_variants;

}

bool kgl::VcfFactory::VcfFileImpl::parseVcfRecord(const std::string& genome_name,
                                                  const seqan::VcfRecord& record,
                                                  ActiveContigMap& active_contig_map,
                                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                  std::shared_ptr<GenomeVariant> genome_variants,
                                                  Phred_t read_quality,
                                                  NucleotideReadCount_t min_read_count,
                                                  double min_proportion) const {

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


//      kgl::Attributes record_attributes;
//  std::cout << "length of record.genotypeInfos: " << seqan::length(record.genotypeInfos) << '\n';
//  std::cout << "record.format: " << toCString(record.format) << '\n';
//  std::cout << "record.ref: " << toCString(record.ref) << '\n';
  std::string info = toCString(record.info);
//  std::cout << "record.info: " << info << '\n';

  // assumes input "key_1=value_1; ...;key_n=value_n"
  std::map<std::string, std::string> info_key_value_map;
  if (not tokenizeVcfInfoKeyValues(info, info_key_value_map)) {

    ExecEnv::log().error("Unable to parse VCF INFO: {}", info);
    return false;

  }

  std::vector<std::pair<char, size_t>> parsed_cigar;
  NucleotideReadCount_t read_depth;
  size_t total_cigar_size;
  std::string cigar;
  double proportion;

  auto result_cigar = info_key_value_map.find(ID_CIGAR_VALUE_);

  if (result_cigar != info_key_value_map.end()) {

    cigar = result_cigar->second;
    if (not parseCigar(cigar, total_cigar_size, parsed_cigar)) {

      return false;

    }

  } else {

    ExecEnv::log().error("VCF factory; cigar field: {} not found in info: {}", ID_CIGAR_VALUE_, info);
    ExecEnv::log().error("VCF file should conform to 'freebayes' format");
    return false;

  }

  auto result_read_depth = info_key_value_map.find(ID_READ_DEPTH_);

  if (result_read_depth != info_key_value_map.end()) {

    read_depth = std::atoll(result_read_depth->second.c_str());

  } else {

    ExecEnv::log().error("VCF factory; read depth field: {} not found in info: {}.", ID_READ_DEPTH_, info);
    ExecEnv::log().error("VCF file should conform to 'freebayes' format");
    return false;

  }

  auto result_proportion = info_key_value_map.find(ID_PROPORTION_);

  if (result_proportion != info_key_value_map.end()) {

    proportion = std::atof(result_read_depth->second.c_str());

  } else {

    ExecEnv::log().error("VCF factory; read depth field: {} not found in info: {}.", ID_PROPORTION_, info);
    ExecEnv::log().error("VCF file should conform to 'freebayes' format");
    return false;

  }

  if (read_depth >= min_read_count and proportion >= min_proportion) {

    std::string reference = seqan::toCString(record.ref);
    std::string alleles = seqan::toCString(record.alt);

    // check sizes.
    if (reference.size() != total_cigar_size) {

      ExecEnv::log().error("VCF factory; reference: {} size: {} does not match cigar: {} size: {}",
                           reference, reference.size(), cigar, total_cigar_size);
      return false;

    }

  }

  return true;

}


bool kgl::VcfFactory::VcfFileImpl::parseCigar(const std::string& cigar,
                                              size_t& total_cigar_size,
                                              std::vector<std::pair<char, size_t>>& parsed_cigar) const {

  parsed_cigar.clear();
  total_cigar_size = 0;
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

    if (it != cigar.end()) {

      switch (*it) {

        case 'X':
        case 'M':
        case 'I':
        case 'D':
          cigar_code = *it;
          break;

        default:
          ExecEnv::log().error("VCF factory; cigar: {} contains unexpected cigar code: {}", cigar, *it);
          return false;

      }

      size_t size = std::atoll(cigar_size.c_str());
      total_cigar_size += size;
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

        std::shared_ptr<ContigFeatures> contig_ptr;
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

    ExecEnv::log().info("This VCF file does not define an allele 'CIGAR' field; this field is required to parse multi-nucleotide alleles");
    ExecEnv::log().info("See the VCF format of the variant caller 'freebayes' for more information.");

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
                                                                  Phred_t read_quality,
                                                                  NucleotideReadCount_t min_read_count,
                                                                  double min_proportion) const {

  return vcf_file_impl_ptr_->readParseVcfFile(genome_name,
                                              genome_db_ptr,
                                              vcf_file_name,
                                              read_quality,
                                              min_read_count,
                                              min_proportion);

}

