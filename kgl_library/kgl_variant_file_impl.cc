//
// Created by kellerberrin on 18/4/20.
//


#include <kel_exec_env.h>
#include "kel_mt_queue.h"
#include "kgl_genome_types.h"
#include "kgl_variant_file.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

#include <string>
#include <vector>
#include <memory>
#include <thread>
#include <fstream>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace kgl = kellerberrin::genome;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::FileVCFIO::commenceIO() {

  raw_io_thread_ptr_ = std::make_unique<std::thread>(&FileVCFIO::rawVCFIO, this);

  for(size_t i = 0; i < PARSER_THREADS_; ++i) {

    vcf_record_thread_vec_.emplace_back(&FileVCFIO::enqueueVCFRecord, this);

  }


}



void kgl::FileVCFIO::rawVCFIO() {

  try {

    std::ifstream vcf_file;

    // Open input file.

    vcf_file.open(fileName());

    if (not vcf_file.good()) {

      ExecEnv::log().critical("I/O error; could not open VCF file: {}", fileName());

    }

    while (true) {

      std::unique_ptr<std::string> line_record_ptr = std::make_unique<std::string>();
      if (std::getline(vcf_file, *line_record_ptr).eof()) {

        // Enqueue the null eof indicator for each consumer thread.
        for(size_t i = 0; i < PARSER_THREADS_; ++i) {

          raw_io_queue_.push(std::unique_ptr<std::string>(nullptr));

        }
        break;

      }

      // Check we have read a string.
      if (line_record_ptr->empty()) {

        ExecEnv::log().error("Empty VCF record string returned");
        continue;

      }

      // Skip header records.
      if ((*line_record_ptr)[0] == HEADER_CHAR_) {

        continue;

      }

      raw_io_queue_.push(std::move(line_record_ptr));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCF file unexpected Seqan I/O exception: {}", e.what());

  }

}


void kgl::FileVCFIO::enqueueVCFRecord() {

  size_t record_line = 0;
  std::unique_ptr<std::string> line_record_ptr;
  while (true) {

    raw_io_queue_.waitAndPop(line_record_ptr);

    if (not line_record_ptr) {

      vcf_record_queue_.push(VcfRecord::EOF_RECORD());
      break;

    }

    ++record_line;

    std::unique_ptr<VcfRecord> vcf_record_ptr(std::make_unique<VcfRecord>());
    if (not parseVCFRecord(line_record_ptr, vcf_record_ptr)) {

      ExecEnv::log().warn("Failed to parse VCF file: {} record line : {}", fileName(), record_line);

    } else {

      vcf_record_queue_.push(std::move(vcf_record_ptr));

    }

  }

}


bool kgl::FileVCFIO::parseVCFRecord( const std::unique_ptr<std::string>& line_record_ptr, const std::unique_ptr<VcfRecord>& vcf_record_ptr) {

  std::vector<std::string> field_vector;
  if (not ParseVCFMiscImpl::tokenize(std::move(*line_record_ptr), VCF_FIELD_DELIMITER_,  field_vector)) {

    ExecEnv::log().error("VCF file: {}, problem parsing record", fileName());
    return false;

  }

  if (field_vector.size() < MINIMUM_VCF_FIELDS_) {

    ExecEnv::log().error("VCF file: {}, record has less than the mandatory field count", fileName());
    return false;

  }

  if (not moveToVcfRecord(field_vector, *vcf_record_ptr)) {

    ExecEnv::log().error("VCF file: {}, cannot parse VCF record field", fileName());
    return false;

  }

  return true;

}


std::unique_ptr<kgl::VcfRecord> kgl::FileVCFIO::readVCFRecord() {

  std::unique_ptr<kgl::VcfRecord> vcf_record;
  vcf_record_queue_.waitAndPop(vcf_record);
  return vcf_record;

}


kgl::VcfHeaderInfo kgl::FileVCFIO::VCFReadHeader() {

  VcfHeaderInfo header_info;

  return header_info;

}


bool kgl::FileVCFIO::moveToVcfRecord(std::vector<std::string>& fields, VcfRecord& vcf_record) {

  try {

      vcf_record.contig_id = std::move(fields[0]);
      vcf_record.offset = std::stoull(fields[1]) - 1; // all offsets are zero based.
      vcf_record.id =  std::move(fields[2]);
      vcf_record.ref =  std::move(fields[3]);
      vcf_record.alt =  std::move(fields[4]);
      vcf_record.qual = std::stod(fields[5]);
      vcf_record.filter = std::move(fields[6]);
      vcf_record.info = std::move(fields[7]);

      if (fields.size() > MINIMUM_VCF_FIELDS_) {

        vcf_record.format = std::move(fields[8]);

        vcf_record.genotypeInfos.clear();
        for (size_t idx = (MINIMUM_VCF_FIELDS_ + 1); idx < fields.size(); ++idx) {

          vcf_record.genotypeInfos.push_back(std::move(fields[idx]));

        }

      }

    return true;

  }
  catch(const std::exception& e) {

    ExecEnv::log().error("Problem parsing record for VCF file: {}, Exception: {} thrown record ignored", e.what(), fileName());
    return false;

  }

}

