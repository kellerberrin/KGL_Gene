//
// Created by kellerberrin on 15/4/20.
//

#include "kel_exec_env.h"
#include "kgl_variant_file.h"

#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

namespace bt = boost;
namespace kgl = kellerberrin::genome;



bool kgl::VCFParseHeader::parseHeader(const std::string& vcf_file_name) {

  std::ifstream vcf_file;

  // Open input file.

  vcf_file.open(vcf_file_name);

  if (not vcf_file.good()) {

    ExecEnv::log().critical("I/O error; could not open VCF file: {}", vcf_file_name);

  }

  try {

    long counter = 0;
    bool found_header = false;


    while (true) {

      std::string record_str;

      if (std::getline(vcf_file, record_str).eof()) break;

      std::string contig_prefix = record_str.substr(0, CONTIG_NAME_FRAGMENT_LENGTH_);
      if (contig_prefix == CONTIG_NAME_FRAGMENT_) {

        std::string contig_string = Utility::trimAllWhiteSpace(record_str);

        ExecEnv::log().info("**VCF Contig Header: {}", record_str);

      }

      std::string line_prefix = record_str.substr(0, FIELD_NAME_FRAGMENT_LENGTH_);
      if (line_prefix == FIELD_NAME_FRAGMENT_) {

        found_header = true;
        size_t field_count = 0;
        bt::char_separator<char> item_key_sep("\t");
        bt::tokenizer<bt::char_separator<char>> tokenize_item(record_str, item_key_sep);
        for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

          if (field_count >= SKIP_FIELD_NAMES_) {

            vcf_genomes_.push_back(*iter_item);

          }

          ++field_count;

        }

        break;

      }

      ++counter;

    }

    vcf_file.close();

    if (not found_header) {

      ExecEnv::log().error("VCF Field Names Not Found");

    } else {

      ExecEnv::log().info("{} Genomes in VCF file {}", vcf_genomes_.size(), vcf_file_name);

    }


  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VVCFParseHeader::parseHeader; VCF file: {}, unexpected I/O exception: {}", vcf_file_name, e.what());

  }

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////



kgl::VcfRecord::VcfRecord(seqan::VcfRecord&& vcf_record, ContigId_t&& contig) : contig_id(contig) {

  offset = vcf_record.beginPos;
  seqan::move(id, vcf_record.id);
  seqan::move(ref, vcf_record.ref);
  seqan::move(alt, vcf_record.alt);
  qual = vcf_record.qual;
  seqan::move(filter, vcf_record.filter);
  seqan::move(info, vcf_record.info);
  seqan::move(format, vcf_record.format);

  auto begin = seqan::begin(vcf_record.genotypeInfos);
  auto end = seqan::end(vcf_record.genotypeInfos);

  for (auto it = begin; it !=end; ++it) {

    std::string genotype;
    seqan::move(genotype, *it);
    genotypeInfos.push_back(std::move(genotype));

  }

}

/*
kgl::VcfRecord::VcfRecord(const seqan::VcfRecord& vcf_record, const ContigId_t& contig) : contig_id(contig) {

  offset = vcf_record.beginPos;
  id = seqan::toCString(vcf_record.id);
  ref = seqan::toCString(vcf_record.ref);
  alt = seqan::toCString(vcf_record.alt);
  qual = vcf_record.qual;
  filter = seqan::toCString(vcf_record.filter);
  info = seqan::toCString(vcf_record.info);
  format = seqan::toCString(vcf_record.format);

  auto begin = seqan::begin(vcf_record.genotypeInfos);
  auto end = seqan::end(vcf_record.genotypeInfos);

  for (auto it = begin; it !=end; ++it) {

    std::string genotype = seqan::toCString(*it);
    genotypeInfos.push_back(std::move(genotype));

  }

}
*/