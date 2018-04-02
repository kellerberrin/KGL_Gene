//
// Created by kellerberrin on 22/01/18.
//

#include "kgl_variant_factory_fbvcf_impl.h"
#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_vcf_cigar_impl.h"
#include "kgl_variant_factory_compound.h"



namespace kgl = kellerberrin::genome;
namespace bt = boost;


bool kgl::FreeBayesVCFImpl::parseVcfRecord(const std::string& genome_name,
                                           const seqan::VcfRecord& record,
                                           std::shared_ptr<const ContigFeatures> contig_ptr,
                                           Phred_t variant_quality,
                                           bool& quality_ok,
                                           size_t& record_variants) {

  std::string info = toCString(record.info);
  // assumes input "key_1=value_1; ...;key_n=value_n"
  std::map<std::string, std::string> info_key_value_map;
  if (not ParseVCFMiscImpl::tokenizeVcfInfoKeyValues(info, info_key_value_map)) {

    ExecEnv::log().error("Unable to parse VCF INFO: {}", info);
    return false;

  }

  std::vector<CigarEditItem> parsed_cigar;
  size_t reference_size;
  size_t alternate_size;
  std::string cigar;
  Phred_t quality = record.qual;

  auto result_cigar = info_key_value_map.find(ID_CIGAR_VALUE_);

  if (result_cigar != info_key_value_map.end()) {

    cigar = result_cigar->second;
    if (not ParseVCFMiscImpl::parseCigar(cigar, reference_size, alternate_size, parsed_cigar)) {

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

    parseCigarItems(genome_name,
                    contig_ptr,
                    parsed_cigar,
                    contig_offset,
                    reference,
                    alternate,
                    quality,
                    info,
                    record_variants);

  } else {

    quality_ok = false;

  }

  return true;

}



