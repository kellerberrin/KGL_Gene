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

    std::optional<std::unique_ptr<VcfRecord>> vcf_record_opt = moveToVcfRecord(std::move(*line_record.value().second));
    if (not vcf_record_opt) {

      ExecEnv::log().warn("FileVCFIO; Failed to parse VCF file: {} record line : {}", fileName(), line_record.value().first);

    } else {

      QueuedVCFRecord queue_record(std::pair<size_t, std::unique_ptr<VcfRecord>>(line_record.value().first, std::move(vcf_record_opt.value())));
      vcf_record_queue_.push(std::move(queue_record));

    }

  }

}


std::optional<std::unique_ptr<VcfRecord>> RecordVCFIO::moveToVcfRecord(std::string&& line_record) {

  try {

    std::unique_ptr<VcfRecord> vcf_record_ptr(std::make_unique<VcfRecord>());

    std::vector<std::string_view> field_views = Utility::view_tokenizer(line_record, VCF_FIELD_DELIMITER_CHAR_);

    if (field_views.size() < MINIMUM_VCF_FIELDS_) {

      ExecEnv::log().error("FileVCFIO; VCF file: {}, record has less than the mandatory field count", fileName());
      return std::nullopt;

    }

    vcf_record_ptr->contig_id = field_views[0];
    vcf_record_ptr->offset = std::stoull(std::string(field_views[1])) - 1; // all offsets are zero based.
    vcf_record_ptr->id = field_views[2];
    vcf_record_ptr->ref = field_views[3];
    // A deletion variant can be signalled by a missing alt value.
    if (field_views[4] == FIELD_NOT_PRESENT_) {

      vcf_record_ptr->alt = "";

    } else {

      vcf_record_ptr->alt = field_views[4];

    }
    // The quality field can be omitted.
    if (field_views[5] == FIELD_NOT_PRESENT_) {

      vcf_record_ptr->qual = 0.0;

    } else {

      vcf_record_ptr->qual = std::stod(std::string(field_views[5]));

    }
    vcf_record_ptr->filter = field_views[6];
    vcf_record_ptr->info = field_views[7];

    if (field_views.size() > MINIMUM_VCF_FIELDS_) {

      vcf_record_ptr->format = field_views[8];

      vcf_record_ptr->genotypeInfos.clear();
      for (size_t idx = (MINIMUM_VCF_FIELDS_ + 1); idx < field_views.size(); ++idx) {

        vcf_record_ptr->genotypeInfos.emplace_back(field_views[idx]);

      }

    }

    // Should not matter since we are moving the string, but just to be on the safe side this goes last.
    vcf_record_ptr->line_record_ptr = std::make_unique<std::string>(std::move(line_record));

    return vcf_record_ptr;

  }
  catch (const std::exception &e) {

    ExecEnv::log().error("FileVCFIO; Problem parsing record for VCF file: {}, Exception: {} thrown; VCF record ignored", fileName(),  e.what());
    ExecEnv::log().error("FileVCFIO; VCF record line: {}", fileName(),  e.what());
    return std::nullopt;

  }

}

} // namespace