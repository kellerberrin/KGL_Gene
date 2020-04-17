//
// Created by kellerberrin on 15/4/20.
//

#include "kel_exec_env.h"
#include "kgl_variant_file.h"

#include <fstream>

#include <boost/tokenizer.hpp>

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

        break; // #CHROM is the last field in the VCF header so stop processing.

      }

      ++counter;

    }

    vcf_file.close();

    if (not found_header) {

      ExecEnv::log().error("VCF Genome Names Not Found");

    } else {

      ExecEnv::log().info("{} Genomes in VCF Header {}, Header lines processed: {}", vcf_genomes_.size(), vcf_file_name, counter);

    }


  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCFParseHeader::parseHeader; VCF file: {}, unexpected I/O exception: {}", vcf_file_name, e.what());

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::ContigId_t kgl::SeqanVCFIO::getContig(int32_t contig_idx) const {

  std::string contig_id = seqan::toCString(contigNames(context(*vcfIn_ptr_))[contig_idx]);

  return contig_id;

}


void kgl::SeqanVCFIO::commenceIO() {

  raw_io_thread_ptr_ = std::make_unique<std::thread>(&SeqanVCFIO::rawVCFIO, this);
  vcf_record_thread_ptr_ = std::make_unique<std::thread>(&SeqanVCFIO::enqueueVCFRecord, this);

}



void kgl::SeqanVCFIO::rawVCFIO() {

  try {

    std::unique_ptr<seqan::VcfRecord> vcf_record_ptr;
    while (true) {

      if (VCFRecordEOF()) {

        raw_io_queue_.push(std::unique_ptr<seqan::VcfRecord>(nullptr));
        break;

      }

      vcf_record_ptr = std::make_unique<seqan::VcfRecord>();
      seqan::readRecord(*vcf_record_ptr, *vcfIn_ptr_);
      raw_io_queue_.push(std::move(vcf_record_ptr));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCF file unexpected Seqan I/O exception: {}", e.what());

  }

}



void kgl::SeqanVCFIO::enqueueVCFRecord() {

  std::unique_ptr<seqan::VcfRecord> vcf_record_ptr;
  while (true) {

    raw_io_queue_.waitAndPop(vcf_record_ptr);

    if (not vcf_record_ptr) {

      vcf_record_queue_.push(VcfRecord::EOF_RECORD());
      break;

    }

    vcf_record_queue_.push(std::make_unique<VcfRecord>(std::move(*vcf_record_ptr), getContig(vcf_record_ptr->rID)));

  }

}


std::unique_ptr<kgl::VcfRecord> kgl::SeqanVCFIO::readVCFRecord() {

  std::unique_ptr<kgl::VcfRecord> vcf_record;
  vcf_record_queue_.waitAndPop(vcf_record);
  return vcf_record;

}


bool kgl::SeqanVCFIO::VCFRecordEOF() { return seqan::atEnd(*vcfIn_ptr_); }

kgl::VcfHeaderInfo kgl::SeqanVCFIO::VCFReadHeader() {

  seqan::VcfHeader vcf_header;
  seqan::readHeader(vcf_header, *vcfIn_ptr_);

  VcfHeaderInfo header_info;
  for (size_t idx = 0; idx != seqan::length(vcf_header); ++idx) {

    std::string key = seqan::toCString(vcf_header[idx].key);
    std::string value = seqan::toCString(vcf_header[idx].value);
    header_info.push_back(std::pair(key,value));

  }

  return header_info;

}

