//
// Created by kellerberrin on 21/4/20.
//



#include <kel_exec_env.h>
#include "kel_utility.h"

#include "kgl_variant_file_vcf_impl.h"

#include <string>
#include <vector>
#include <memory>

// Implementation file classes need to be defined within namespaces.
namespace kellerberrin::genome {



//////////////////////////////////////////////////////////////////////////////////////////////////////////////

RecordVCFIO::~RecordVCFIO() noexcept {

  // Ensure all worker threads are shutdown, suppress exceptions.
  try {

    launch_token_.wait();

  } catch(...) {}

}

bool RecordVCFIO::commenceVCFIO(const std::string& vcf_file_name) {

  if (not file_data_.commenceIO(vcf_file_name)) {

    // Enqueue an eof marker further up the pipeline.
    enqueueEOF();
    return false;

  }

  launch_token_ = detached_launch_.enqueueTask(&RecordVCFIO::launchThreads, this);

  return true;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Launch the VCF record parser worker threads and await completion.

void RecordVCFIO::launchThreads() {

  ThreadPool vcf_record_threads{PARSER_THREADS_};
  std::vector<std::future<void>> thread_futures;
  for (size_t index = 0; index < vcf_record_threads.threadCount(); ++index) {

    thread_futures.push_back(vcf_record_threads.enqueueTask(&RecordVCFIO::enqueueVCFRecord, this));

  }

  // Wait until processing is complete.
  for (auto const& future : thread_futures) {

    future.wait();

  }

  // Enqueue an eof marker further up the pipeline.
  enqueueEOF();

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parse the line record into a VCF record and enqueue the record.

void RecordVCFIO::enqueueVCFRecord() {

  while (true) {

    IOLineRecord line_record = file_data_.readIORecord();
    if (not line_record) { // check for EOF condition.

      // push the eof marker back on the queue.
      file_data_.enqueueEOF();
      break;

    }

    // Skip header records by examining the first char for '#'
    std::string& line_string = *line_record.value().second;
    if (line_string[0] == HEADER_CHAR_) {

      continue;

    }

    std::optional<std::unique_ptr<VcfRecord>> vcf_record_opt = moveToVcfRecord(std::move(*line_record.value().second));
    if (not vcf_record_opt) {

      ExecEnv::log().warn("FileVCFIO; Failed to parse VCF file: {} record line : {}", file_data_.fileName(), line_record.value().first);

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

      ExecEnv::log().error("FileVCFIO; VCF file: {}, record has less than the mandatory field count", file_data_.fileName());
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

    ExecEnv::log().error("FileVCFIO; Problem parsing record for VCF file: {}, Exception: {} thrown; VCF record ignored", file_data_.fileName(),  e.what());
    ExecEnv::log().error("FileVCFIO; VCF record line: {}", file_data_.fileName(),  e.what());
    return std::nullopt;

  }

}

} // namespace