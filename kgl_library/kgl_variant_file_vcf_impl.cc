//
// Created by kellerberrin on 21/4/20.
//



#include <kel_exec_env.h>

#include "kgl_variant_file_vcf_impl.h"

#include <string>
#include <vector>
#include <memory>
#include <thread>

// Implementation file classes need to be defined within namespaces.
namespace kellerberrin::genome {



//////////////////////////////////////////////////////////////////////////////////////////////////////////////

RecordVCFIO::~RecordVCFIO() {

  for (auto &thread : vcf_record_thread_vec_) {

    thread.join();

  }

}

void RecordVCFIO::commenceVCFIO(size_t reader_threads) {

  reader_threads_ = reader_threads;

  commenceIO(PARSER_THREADS_);

  for (size_t i = 0; i < PARSER_THREADS_; ++i) {

    vcf_record_thread_vec_.emplace_back(&RecordVCFIO::enqueueVCFRecord, this);

  }



}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pass in the text or gzip stream as an rvalue to a pointer.

void RecordVCFIO::enqueueVCFRecord() {

  while (true) {

    IOLineRecord line_record = readIORecord();

    if (not line_record) { // check for EOF condition.

      // Only queue as many EOF (null optional) tokens as the number of reader threads.
      for (size_t i = 0; i < reader_threads_; ++i) {

        vcf_record_queue_.push(std::nullopt);

      }
      break;

    }

    // Skip header records by examining the first char for '#'
    std::string& line_string = *line_record.value().second;
    if (line_string[0] == HEADER_CHAR_) {

      continue;

    }

    std::unique_ptr<VcfRecord> vcf_record_ptr(std::make_unique<VcfRecord>());
    if (not parseVCFRecord(std::move(line_record.value().second), vcf_record_ptr)) {

      ExecEnv::log().warn("FileVCFIO; Failed to parse VCF file: {} record line : {}", fileName(), line_record.value().first);

    } else {

      QueuedVCFRecord queue_record(std::pair<size_t, std::unique_ptr<VcfRecord>>(line_record.value().first, std::move(vcf_record_ptr)));
      vcf_record_queue_.push(std::move(queue_record));

    }

  }

}


bool RecordVCFIO::parseVCFRecord( std::unique_ptr<const std::string> line_record_ptr,
                                  const std::unique_ptr<VcfRecord> &vcf_record_ptr) {

  std::vector<std::string> field_vector = Utility::char_tokenizer(*line_record_ptr, VCF_FIELD_DELIMITER_CHAR_);

  if (field_vector.size() < MINIMUM_VCF_FIELDS_) {

    ExecEnv::log().error("FileVCFIO; VCF file: {}, record has less than the mandatory field count", fileName());
    return false;

  }

  if (not moveToVcfRecord(std::move(line_record_ptr),field_vector, *vcf_record_ptr)) {

    ExecEnv::log().error("FileVCFIO; VCF file: {}, cannot parse VCF record field", fileName());
    return false;

  }

  return true;

}



bool RecordVCFIO::moveToVcfRecord(std::unique_ptr<const std::string> line_record_ptr, std::vector<std::string> &fields, VcfRecord &vcf_record) {

  try {

    vcf_record.line_record_ptr = std::move(line_record_ptr);

    vcf_record.contig_id = fields[0];
    vcf_record.offset = std::stoull(fields[1]) - 1; // all offsets are zero based.
    vcf_record.id = std::move(fields[2]);
    vcf_record.ref = std::move(fields[3]);
    // A deletion variant can be signalled by a missing alt value.
    if (fields[4] == FIELD_NOT_PRESENT_) {

      vcf_record.alt = "";

    } else {

      vcf_record.alt = std::move(fields[4]);

    }
    // The quality field can be omitted.
    if (fields[5] == FIELD_NOT_PRESENT_) {

      vcf_record.qual = 0.0;

    } else {

      vcf_record.qual = std::stod(fields[5]);

    }
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
  catch (const std::exception &e) {

    ExecEnv::log().error("FileVCFIO; Problem parsing record for VCF file: {}, Exception: {} thrown; VCF record ignored", fileName(),  e.what());
    ExecEnv::log().error("FileVCFIO; VCF record line: {}", fileName(),  e.what());
    return false;

  }

}

} // namespace